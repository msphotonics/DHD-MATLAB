"""Replace image-filled tiling patterns (as written by MATLAB exportgraphics
'vector') with direct image draws, so viewers such as macOS Preview, which
mis-phase tiling patterns, render the images correctly.
Usage: python fix_pdf_patterns.py in.pdf out.pdf   (requires: pip install pypdf)
"""
import sys
from pypdf import PdfReader, PdfWriter
from pypdf.generic import (ContentStream, NameObject, NumberObject, FloatObject,
                           DictionaryObject, ArrayObject)

def mul(a, b):  # PDF matrices [a b c d e f]; returns a x b
    return [a[0]*b[0]+a[1]*b[2], a[0]*b[1]+a[1]*b[3],
            a[2]*b[0]+a[3]*b[2], a[2]*b[1]+a[3]*b[3],
            a[4]*b[0]+a[5]*b[2]+b[4], a[4]*b[1]+a[5]*b[3]+b[5]]
def inv(m):
    a,b,c,d,e,f = m; det = a*d-b*c
    ia,ib,ic,id_ = d/det, -b/det, -c/det, a/det
    return [ia, ib, ic, id_, -(e*ia+f*ic), -(e*ib+f*id_)]
I = [1,0,0,1,0,0]
PATH = {b'm', b'l', b'c', b'v', b'y', b'h', b're'}

def image_pattern(pat):
    """Return (image name, image ref, matrix of image in pattern space) or None."""
    pat = pat.get_object()
    if int(pat.get('/PatternType', 0)) != 1: return None
    xo = pat.get('/Resources', {}).get('/XObject', {})
    ops = ContentStream(pat, None).operations
    ctm, name = I, None
    for args, op in ops:
        if op == b'cm': ctm = mul([float(x) for x in args], ctm)
        elif op == b'Do': name = args[0]
        elif op in (b'q', b'Q'): pass
        else: return None
    if name is None or name not in xo: return None
    bbox = [float(x) for x in pat['/BBox']]
    if bbox[0] != 0 or bbox[1] != 0 or abs(float(pat['/XStep'])-bbox[2]) > 1e-6: return None
    return name, xo.raw_get(name), ctm, [float(x) for x in pat.get('/Matrix', I)]

def fix_page(page, writer):
    res = page['/Resources']
    pats = res.get('/Pattern', {})
    info = {}
    for n in pats:
        r = image_pattern(pats[n])
        if r: info[n] = r
    if not info: return 0
    cs = ContentStream(page.get_contents(), page.pdf)
    out, stack, ctm, cur, n_fix = [], [], I, None, 0
    xobj = res.get('/XObject')
    if xobj is None:
        xobj = DictionaryObject(); res[NameObject('/XObject')] = xobj
    else:
        xobj = xobj.get_object()
    for args, op in cs.operations:
        if op == b'q': stack.append((ctm, cur))
        elif op == b'Q': ctm, cur = stack.pop() if stack else (I, None)
        elif op == b'cm': ctm = mul([float(x) for x in args], ctm)
        elif op == b'scn' and args and isinstance(args[-1], NameObject):
            cur = args[-1] if args[-1] in info else None
        elif op in (b'f', b'f*', b'F') and cur is not None:
            name, ref, imgm, patm = info[cur]
            newname = NameObject('/UP' + cur[1:] )
            xobj[newname] = ref
            # image space -> pattern space -> default space -> current user space
            m = mul(mul(imgm, patm), inv(ctm))
            out.append(([], b'W*' if op == b'f*' else b'W'))
            out.append(([], b'n'))
            out.append(([], b'q'))
            out.append(([FloatObject(round(x, 6)) for x in m], b'cm'))
            out.append(([newname], b'Do'))
            out.append(([], b'Q'))
            n_fix += 1
            continue
        out.append((args, op))
    cs.operations = out
    page.replace_contents(cs)
    return n_fix

r = PdfReader(sys.argv[1]); w = PdfWriter(clone_from=r)
n = sum(fix_page(p, w) for p in w.pages)
w.write(sys.argv[2]); print(f'{sys.argv[1]}: replaced {n} pattern fills')
