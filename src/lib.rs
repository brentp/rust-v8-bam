use std::ffi::c_void;
use std::sync::Once;

use anyhow::{Result, anyhow};
use rust_htslib::bam;
use rust_htslib::bam::record::{Aux, Cigar};

use v8::{self, Global};

static INIT_V8: Once = Once::new();

fn init_v8_once() {
    INIT_V8.call_once(|| {
        let platform = v8::new_default_platform(0, false).make_shared();
        v8::V8::initialize_platform(platform);
        v8::V8::initialize();
    });
}

/// Engine that owns a V8 isolate, context, compiled filter function,
/// and a reusable `aln` object.
pub struct JsBamFilterEngine {
    isolate: v8::OwnedIsolate,
    context: Global<v8::Context>,
    filter_fn: Global<v8::Function>,
    aln_obj: Global<v8::Object>,
}

impl JsBamFilterEngine {
    /// Create a new engine with a JS filter expression or body.
    ///
    /// `expr` can be:
    ///   "aln.mapq > 10 && aln.qname.startsWith('q23')"
    /// or:
    ///   "return aln.mapq > 10 && aln.qname.startsWith('q23');"
    pub fn new(expr: &str) -> Result<Self> {
        init_v8_once();

        let mut isolate = v8::Isolate::new(Default::default());

        // Create locals first, then convert to globals
        let (ctx_global, filter_global, aln_obj_global) = {
            // Pinned handle scope
            v8::scope!(let hs, &mut isolate);

            // Context
            let context = v8::Context::new(hs, Default::default());
            v8::scope_with_context!(let scope, hs, context);

            // Build full JS source: define `filter(aln)` and helper function(s)
            let source = make_filter_source(expr);
            let filter_fn = compile_filter_function(scope, context, &source)?;

            // Make aln template (lazy accessors mapq, qname, flag, pos)
            let aln_tmpl = make_aln_template(scope);

            // Create a single reusable aln object from the template
            let aln_obj = aln_tmpl
                .new_instance(scope)
                .ok_or_else(|| anyhow!("failed to create aln object"))?;

            // Install global Rust helpers into the context (e.g. hasFlag)
            install_rust_helpers(scope, context);

            // Convert to globals
            let ctx_global = Global::new(scope, context);
            let filter_global = Global::new(scope, filter_fn);
            let aln_obj_global = Global::new(scope, aln_obj);

            (ctx_global, filter_global, aln_obj_global)
        };

        Ok(Self {
            isolate,
            context: ctx_global,
            filter_fn: filter_global,
            aln_obj: aln_obj_global,
        })
    }

    /// Run the JS filter on a single BAM record.
    pub fn record_passes(&mut self, rec: &bam::Record, header: &bam::HeaderView) -> Result<bool> {
        v8::scope!(let hs, &mut self.isolate);
        let context = v8::Local::new(hs, &self.context);
        v8::scope_with_context!(let scope, hs, context);

        let filter_fn = v8::Local::new(scope, &self.filter_fn);
        let aln_obj = v8::Local::new(scope, &self.aln_obj);

        let hdr_ptr = header as *const bam::HeaderView as *mut c_void;
        aln_obj.set_aligned_pointer_in_internal_field(0, hdr_ptr);

        // Store pointer to `bam::Record` in internal field 0.
        // Lifetime: only valid during this call.
        let ptr = rec as *const bam::Record as *mut c_void;
        aln_obj.set_aligned_pointer_in_internal_field(1, ptr);

        let undefined = v8::undefined(scope).into();
        let args = [aln_obj.into()];
        let result = filter_fn
            .call(scope, undefined, &args)
            .ok_or_else(|| anyhow!("filter() threw or returned empty"))?;

        Ok(result.boolean_value(scope))
    }
}

/// Build the JS source that defines the filter.
fn make_filter_source(user_expr: &str) -> String {
    // Allow "and"/"or" as sugar
    let mut expr = user_expr.to_string();
    expr = expr.replace(" and ", " && ");
    expr = expr.replace(" or ", " || ");

    let body = if expr.contains("return") {
        expr
    } else {
        format!("return {};", expr)
    };

    // We also expose Rust helpers (installed separately as globals),
    // e.g. hasFlag(flag, mask).
    //
    // This string only needs to define filter(aln).
    format!(
        r#"
        function filter(aln) {{
            {body}
        }}
        "#,
        body = body
    )
}

/// Compile `filter(aln)` and return the function handle.
fn compile_filter_function<'s>(
    scope: &mut v8::ContextScope<'s, '_, v8::HandleScope<'_>>,
    context: v8::Local<'s, v8::Context>,
    source: &str,
) -> Result<v8::Local<'s, v8::Function>> {
    let code = v8::String::new(scope, source)
        .ok_or_else(|| anyhow!("failed to create JS source string"))?;
    let script =
        v8::Script::compile(scope, code, None).ok_or_else(|| anyhow!("failed to compile JS"))?;
    script
        .run(scope)
        .ok_or_else(|| anyhow!("failed to run JS"))?;

    let global = context.global(scope);
    let name = v8::String::new(scope, "filter").unwrap().into();
    let value = global
        .get(scope, name)
        .ok_or_else(|| anyhow!("global.filter not found"))?;
    let func = v8::Local::<v8::Function>::try_from(value)
        .map_err(|_| anyhow!("filter is not a function"))?;
    Ok(func)
}

/// Create an ObjectTemplate for `aln` with lazy accessors:
/// for chrom, mapq, qname, flag, pos, start, end, aux(tag), etc
fn make_aln_template<'s>(
    scope: &mut v8::ContextScope<'s, '_, v8::HandleScope<'_>>,
) -> v8::Local<'s, v8::ObjectTemplate> {
    let tmpl = v8::ObjectTemplate::new(scope);
    // 0: header view, 1: record
    tmpl.set_internal_field_count(2);

    let mapq = v8::String::new(scope, "mapq").unwrap();
    tmpl.set_accessor(mapq.into(), aln_mapq_getter);

    let qname = v8::String::new(scope, "qname").unwrap();
    tmpl.set_accessor(qname.into(), aln_qname_getter);

    let flag = v8::String::new(scope, "flag").unwrap();
    tmpl.set_accessor(flag.into(), aln_flag_getter);

    let pos = v8::String::new(scope, "pos").unwrap();
    tmpl.set_accessor(pos.into(), aln_pos_getter);

    let start = v8::String::new(scope, "start").unwrap();
    tmpl.set_accessor(start.into(), aln_pos_getter);
    let end = v8::String::new(scope, "end").unwrap();
    tmpl.set_accessor(end.into(), aln_end_getter);

    let chrom = v8::String::new(scope, "chrom").unwrap();
    tmpl.set_accessor(chrom.into(), aln_chrom_getter);

    let cigar = v8::String::new(scope, "cigar").unwrap();
    tmpl.set_accessor(cigar.into(), aln_cigar_getter);

    // Add aux(tag) method
    let aux_fn = v8::FunctionTemplate::new(scope, aln_aux_method);
    let aux_name = v8::String::new(scope, "aux").unwrap();
    tmpl.set(aux_name.into(), aux_fn.into());

    tmpl
}

/// Install global helper functions implemented in Rust.
/// Example: `hasFlag(flag, mask)` → boolean.
fn install_rust_helpers(
    scope: &mut v8::ContextScope<'_, '_, v8::HandleScope<'_>>,
    context: v8::Local<v8::Context>,
) {
    let global = context.global(scope);

    // hasFlag(flag, mask) => (flag & mask) != 0
    let name = v8::String::new(scope, "hasFlag").unwrap();
    let func = v8::Function::new(scope, has_flag_callback).unwrap();
    global.set(scope, name.into(), func.into());
}

#[inline(always)]
fn record_from_obj<'s>(obj: v8::Local<v8::Object>) -> &'s bam::Record {
    let ptr = unsafe { obj.get_aligned_pointer_from_internal_field(1) } as *const bam::Record;
    unsafe { &*ptr }
}

#[inline(always)]
fn header_from_obj<'s>(obj: v8::Local<v8::Object>) -> &'s bam::HeaderView {
    let ptr = unsafe { obj.get_aligned_pointer_from_internal_field(0) } as *const bam::HeaderView;
    unsafe { &*ptr as &bam::HeaderView }
}

// ========== Accessors: aln.mapq, aln.qname, aln.flag, aln.pos ==========

#[allow(clippy::needless_pass_by_value)]
fn aln_mapq_getter(
    scope: &mut v8::PinScope,
    _name: v8::Local<v8::Name>,
    args: v8::PropertyCallbackArguments,
    mut rv: v8::ReturnValue,
) {
    let this = args.this();
    let rec = record_from_obj(this);
    let v = v8::Integer::new_from_unsigned(scope, rec.mapq() as u32);
    rv.set(v.into());
}

#[allow(clippy::needless_pass_by_value)]
fn aln_qname_getter(
    scope: &mut v8::PinScope,
    _name: v8::Local<v8::Name>,
    args: v8::PropertyCallbackArguments,
    mut rv: v8::ReturnValue,
) {
    let this = args.this();
    let rec = record_from_obj(this);
    let qname_bytes = rec.qname();
    let qname = std::str::from_utf8(qname_bytes).unwrap_or("");
    let s = v8::String::new(scope, qname).unwrap();
    rv.set(s.into());
}

#[allow(clippy::needless_pass_by_value)]
fn aln_flag_getter(
    scope: &mut v8::PinScope,
    _name: v8::Local<v8::Name>,
    args: v8::PropertyCallbackArguments,
    mut rv: v8::ReturnValue,
) {
    let this = args.this();
    let rec = record_from_obj(this);
    let flags = rec.flags() as u32;
    let v = v8::Integer::new_from_unsigned(scope, flags);
    rv.set(v.into());
}

fn aln_chrom_getter(
    scope: &mut v8::PinScope,
    _name: v8::Local<v8::Name>,
    args: v8::PropertyCallbackArguments,
    mut rv: v8::ReturnValue,
) {
    let this = args.this();
    let rec = record_from_obj(this);
    let header = header_from_obj(this);
    let tid = rec.tid() as u32;
    let chrom = header.tid2name(tid);
    let chrom = std::str::from_utf8(chrom).unwrap_or("");
    let s = v8::String::new(scope, chrom).unwrap();
    rv.set(s.into());
}

fn aln_end_getter(
    scope: &mut v8::PinScope,
    _name: v8::Local<v8::Name>,
    args: v8::PropertyCallbackArguments,
    mut rv: v8::ReturnValue,
) {
    let this = args.this();
    let rec = record_from_obj(this);
    let end = rec.cigar().end_pos() as u32;
    let v = v8::Integer::new_from_unsigned(scope, end);
    rv.set(v.into());
}

#[allow(clippy::needless_pass_by_value)]
fn aln_pos_getter(
    scope: &mut v8::PinScope,
    _name: v8::Local<v8::Name>,
    args: v8::PropertyCallbackArguments,
    mut rv: v8::ReturnValue,
) {
    let this = args.this();
    let rec = record_from_obj(this);
    // BAM pos is 0-based; we expose that directly.
    let pos = rec.pos() as i32;
    let v = v8::Integer::new(scope, pos);
    rv.set(v.into());
}

fn aln_cigar_getter(
    scope: &mut v8::PinScope,
    _name: v8::Local<v8::Name>,
    args: v8::PropertyCallbackArguments,
    mut rv: v8::ReturnValue,
) {
    let this = args.this();
    let rec = record_from_obj(this);
    let cigar = rec.cigar();
    let js_arr = v8::Array::new(scope, cigar.len() as i32);
    let len_key = v8::String::new(scope, "length").unwrap();
    let op_key = v8::String::new(scope, "op").unwrap();
    let cref_key = v8::String::new(scope, "consumes_ref").unwrap();
    let cquery_key = v8::String::new(scope, "consumes_query").unwrap();

    for (i, op) in cigar.iter().enumerate() {
        let (op_name, consumes_ref, consumes_query, len) = cigar_op_info(op);

        let obj = v8::Object::new(scope);

        let len_val = v8::Integer::new_from_unsigned(scope, len);
        let op_val = v8::String::new(scope, op_name).unwrap();
        let cref_val = v8::Boolean::new(scope, consumes_ref);
        let cquery_val = v8::Boolean::new(scope, consumes_query);

        obj.set(scope, len_key.into(), len_val.into());
        obj.set(scope, op_key.into(), op_val.into());
        obj.set(scope, cref_key.into(), cref_val.into());
        obj.set(scope, cquery_key.into(), cquery_val.into());

        js_arr.set_index(scope, i as u32, obj.into());
    }

    rv.set(js_arr.into());
}

fn cigar_op_info(op: &Cigar) -> (&'static str, bool, bool, u32) {
    match *op {
        Cigar::Match(len) => ("Match", true, true, len),
        Cigar::Ins(len) => ("Ins", false, true, len),
        Cigar::Del(len) => ("Del", true, false, len),
        Cigar::RefSkip(len) => ("RefSkip", true, false, len),
        Cigar::SoftClip(len) => ("SoftClip", false, true, len),
        Cigar::HardClip(len) => ("HardClip", false, false, len),
        Cigar::Pad(len) => ("Pad", false, false, len),
        Cigar::Equal(len) => ("Equal", true, true, len),
        Cigar::Diff(len) => ("Diff", true, true, len),
    }
}

// ========== Method: aln.aux(tag) ==========

#[allow(clippy::needless_pass_by_value)]
fn aln_aux_method(
    scope: &mut v8::PinScope,
    args: v8::FunctionCallbackArguments,
    mut rv: v8::ReturnValue,
) {
    // Get the aln object (this) and extract the BAM record
    let this = args.this();
    let rec = record_from_obj(this);

    // Get tag name from first argument
    let tag_arg = args.get(0);
    if !tag_arg.is_string() {
        rv.set(v8::null(scope).into());
        return;
    }
    let tag_str = tag_arg.to_rust_string_lossy(scope);
    let tag_bytes = tag_str.as_bytes();

    // Tag must be exactly 2 characters
    if tag_bytes.len() != 2 {
        rv.set(v8::null(scope).into());
        return;
    }

    // Look up the aux tag
    match rec.aux(tag_bytes) {
        Ok(aux) => {
            let js_val = aux_to_js_value(scope, aux);
            rv.set(js_val);
        }
        Err(_) => {
            // Tag not found - return null
            rv.set(v8::null(scope).into());
        }
    }
}

/// Convert a rust_htslib Aux value to a V8 value
fn aux_to_js_value<'s, 'i>(
    scope: &mut v8::PinScope<'s, 'i>,
    aux: Aux<'_>,
) -> v8::Local<'s, v8::Value> {
    match aux {
        Aux::I8(v) => v8::Integer::new(scope, v as i32).into(),
        Aux::U8(v) => v8::Integer::new_from_unsigned(scope, v as u32).into(),
        Aux::I16(v) => v8::Integer::new(scope, v as i32).into(),
        Aux::U16(v) => v8::Integer::new_from_unsigned(scope, v as u32).into(),
        Aux::I32(v) => v8::Integer::new(scope, v).into(),
        Aux::U32(v) => v8::Integer::new_from_unsigned(scope, v).into(),
        Aux::Float(v) => v8::Number::new(scope, v as f64).into(),
        Aux::Double(v) => v8::Number::new(scope, v).into(),
        Aux::Char(v) => {
            let bytes = [v];
            let s = String::from_utf8_lossy(&bytes);
            v8::String::new(scope, &s).unwrap().into()
        }
        Aux::String(s) => v8::String::new(scope, s).unwrap().into(),
        Aux::HexByteArray(hex_view) => {
            // HexByteArray is already a hex-encoded string view
            v8::String::new(scope, hex_view.as_ref()).unwrap().into()
        }
        // Array types - return as JS arrays
        Aux::ArrayI8(arr) => {
            let js_arr = v8::Array::new(scope, arr.iter().count() as i32);
            for (i, v) in arr.iter().enumerate() {
                let js_v = v8::Integer::new(scope, v as i32);
                js_arr.set_index(scope, i as u32, js_v.into());
            }
            js_arr.into()
        }
        Aux::ArrayU8(arr) => {
            let js_arr = v8::Array::new(scope, arr.iter().count() as i32);
            for (i, v) in arr.iter().enumerate() {
                let js_v = v8::Integer::new_from_unsigned(scope, v as u32);
                js_arr.set_index(scope, i as u32, js_v.into());
            }
            js_arr.into()
        }
        Aux::ArrayI16(arr) => {
            let js_arr = v8::Array::new(scope, arr.iter().count() as i32);
            for (i, v) in arr.iter().enumerate() {
                let js_v = v8::Integer::new(scope, v as i32);
                js_arr.set_index(scope, i as u32, js_v.into());
            }
            js_arr.into()
        }
        Aux::ArrayU16(arr) => {
            let js_arr = v8::Array::new(scope, arr.iter().count() as i32);
            for (i, v) in arr.iter().enumerate() {
                let js_v = v8::Integer::new_from_unsigned(scope, v as u32);
                js_arr.set_index(scope, i as u32, js_v.into());
            }
            js_arr.into()
        }
        Aux::ArrayI32(arr) => {
            let js_arr = v8::Array::new(scope, arr.iter().count() as i32);
            for (i, v) in arr.iter().enumerate() {
                let js_v = v8::Integer::new(scope, v);
                js_arr.set_index(scope, i as u32, js_v.into());
            }
            js_arr.into()
        }
        Aux::ArrayU32(arr) => {
            let js_arr = v8::Array::new(scope, arr.iter().count() as i32);
            for (i, v) in arr.iter().enumerate() {
                let js_v = v8::Integer::new_from_unsigned(scope, v);
                js_arr.set_index(scope, i as u32, js_v.into());
            }
            js_arr.into()
        }
        Aux::ArrayFloat(arr) => {
            let js_arr = v8::Array::new(scope, arr.iter().count() as i32);
            for (i, v) in arr.iter().enumerate() {
                let js_v = v8::Number::new(scope, v as f64);
                js_arr.set_index(scope, i as u32, js_v.into());
            }
            js_arr.into()
        }
    }
}

// ========== Rust helper: hasFlag(flag, mask) ==========

#[allow(clippy::needless_pass_by_value)]
fn has_flag_callback(
    scope: &mut v8::PinScope,
    args: v8::FunctionCallbackArguments,
    mut rv: v8::ReturnValue,
) {
    let flag_val = args.get(0);
    let mask_val = args.get(1);

    let flag = flag_val.integer_value(scope).unwrap_or(0) as u32;
    let mask = mask_val.integer_value(scope).unwrap_or(0) as u32;

    let result = (flag & mask) != 0;
    let js_bool = v8::Boolean::new(scope, result);
    rv.set(js_bool.into());
}
