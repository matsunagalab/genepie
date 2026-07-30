"""Small DCD rewriting helpers for lazy-reader regression tests."""

from pathlib import Path
import struct


def rewrite_dcd(
    source,
    destination,
    *,
    output_endian: str = ">",
    box_values=None,
):
    """Rewrite a standard fixed-atom DCD with another endian or unit cell."""
    data = Path(source).read_bytes()
    if struct.unpack_from("<i", data, 0)[0] == 84:
        input_endian = "<"
    elif struct.unpack_from(">i", data, 0)[0] == 84:
        input_endian = ">"
    else:
        raise ValueError("Not a supported DCD file")

    offset = 0

    def read_record():
        nonlocal offset
        size = struct.unpack_from(f"{input_endian}i", data, offset)[0]
        offset += 4
        payload = data[offset:offset + size]
        offset += size
        end_size = struct.unpack_from(f"{input_endian}i", data, offset)[0]
        offset += 4
        if end_size != size:
            raise ValueError("Invalid DCD record")
        return payload

    def pack_record(payload):
        size = len(payload)
        return (
            struct.pack(f"{output_endian}i", size)
            + payload
            + struct.pack(f"{output_endian}i", size)
        )

    header = read_record()
    signature = header[:4]
    controls = struct.unpack(f"{input_endian}20i", header[4:])
    output = [
        pack_record(signature + struct.pack(f"{output_endian}20i", *controls))
    ]

    title = read_record()
    ntitle = struct.unpack(f"{input_endian}i", title[:4])[0]
    output.append(
        pack_record(
            struct.pack(f"{output_endian}i", ntitle) + title[4:]
        )
    )

    atom_record = read_record()
    natom = struct.unpack(f"{input_endian}i", atom_record)[0]
    output.append(pack_record(struct.pack(f"{output_endian}i", natom)))

    nframe = controls[0]
    coordinate_bytes = 3 * (8 + 4 * natom)
    remaining = len(data) - offset
    has_box = remaining == nframe * (coordinate_bytes + 56)
    if not has_box and remaining != nframe * coordinate_bytes:
        raise ValueError("Unsupported DCD frame layout")

    for _ in range(nframe):
        if has_box:
            box = struct.unpack(f"{input_endian}6d", read_record())
            if box_values is not None:
                box = tuple(box_values)
            output.append(
                pack_record(struct.pack(f"{output_endian}6d", *box))
            )
        for _axis in range(3):
            values = struct.unpack(f"{input_endian}{natom}f", read_record())
            output.append(
                pack_record(
                    struct.pack(f"{output_endian}{natom}f", *values)
                )
            )

    Path(destination).write_bytes(b"".join(output))
    return nframe, natom, has_box

