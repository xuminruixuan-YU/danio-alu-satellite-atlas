"""Repair omitted PNG predictor bit depth without changing image bytes or vectors."""
from pathlib import Path
import hashlib
import zlib
from pypdf import PdfReader, PdfWriter
from pypdf.generic import DictionaryObject, NameObject, NumberObject


def image_objects(reader):
    """Yield images and soft masks, including nested Form resources, once."""
    seen = set()
    def walk(resources, page, prefix=''):
        for name, ref in resources.get('/XObject', {}).items():
            obj = ref.get_object()
            if id(obj) in seen:
                continue
            seen.add(id(obj))
            if obj.get('/Subtype') == '/Form':
                yield from walk(obj.get('/Resources', {}), page, prefix + str(name) + '/')
            elif obj.get('/Subtype') == '/Image':
                yield page, prefix + str(name), obj
                mask = obj.get('/SMask')
                if mask is not None:
                    mask = mask.get_object()
                    if id(mask) not in seen:
                        seen.add(id(mask))
                        yield page, prefix + str(name) + ':SMask', mask
    for number, page in enumerate(reader.pages, 1):
        yield from walk(page.get('/Resources', {}), number)


def inspect_predictors(reader):
    problems = []
    for page, name, im in image_objects(reader):
        params = im.get('/DecodeParms')
        if im.get('/Filter') != '/FlateDecode' or params is None:
            continue
        params = params.get_object()
        if int(params.get('/Predictor', 1)) < 10:
            continue
        bits = int(im['/BitsPerComponent'])
        declared = int(params.get('/BitsPerComponent', 8))
        columns = int(params.get('/Columns', 1))
        colors = int(params.get('/Colors', 1))
        height = int(im['/Height'])
        raw = zlib.decompress(im._data)
        expected = ((columns * colors * declared + 7) // 8 + 1) * height
        correct = ((columns * colors * bits + 7) // 8 + 1) * height
        if len(raw) != expected or declared != bits:
            if columns != int(im['/Width']) or len(raw) != correct:
                raise ValueError(f'Unexpected image encoding at page {page} {name}; manual diagnosis needed')
            row_bytes = (columns * colors * bits + 7) // 8
            assert all(raw[i] <= 4 for i in range(0, len(raw), row_bytes + 1))
            problems.append((im, dict(page=page, image=name, image_bits=bits,
                        declared_predictor_bits=declared, raw_bytes=len(raw))))
    return problems


def repair_pdf(path, destination=None):
    path = Path(path)
    reader = PdfReader(path)
    problems = inspect_predictors(reader)
    if not problems:
        return []
    before_contents = [hashlib.sha256(p.get_contents().get_data()).hexdigest() for p in reader.pages]
    before_images = {(p,n):hashlib.sha256(im._data).hexdigest() for p,n,im in image_objects(reader)}
    for im, info in problems:
        params = DictionaryObject(im['/DecodeParms'].get_object())
        params[NameObject('/BitsPerComponent')] = NumberObject(info['image_bits'])
        im[NameObject('/DecodeParms')] = params
    writer = PdfWriter(clone_from=reader)
    destination = Path(destination) if destination else path
    temporary = destination.with_suffix('.repairing.pdf')
    with temporary.open('wb') as stream:
        writer.write(stream)
    result = PdfReader(temporary)
    assert not inspect_predictors(result)
    assert before_contents == [hashlib.sha256(p.get_contents().get_data()).hexdigest() for p in result.pages]
    assert before_images == {(p,n):hashlib.sha256(im._data).hexdigest() for p,n,im in image_objects(result)}
    # Require actual decoder success as well as plausible compressed lengths.
    for _, _, image in image_objects(result):
        data = image.get_data()
        cs = image.get('/ColorSpace')
        components = 3 if cs == '/DeviceRGB' else 4 if cs == '/DeviceCMYK' else 1
        expected = ((int(image['/Width']) * components * int(image['/BitsPerComponent']) + 7) // 8) * int(image['/Height'])
        assert len(data) == expected, (len(data), expected)
    temporary.replace(destination)
    return [info for _,info in problems]
