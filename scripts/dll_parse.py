import os


parse_command = {
    # My dumpbin dir
    "win": "(\"D:\\Apps\\Microsoft Visual Studio\\2019\\Community\\VC\\Tools\\MSVC\\14.29.30133\\bin\\Hostx64\\x64"
           "\\dumpbin\" -EXPORTS {dll}.dll)",
    "linux": "nm -D {dll}.so"
}


def parse(dll):
    if os.name == "nt":
        # Windows
        os.system(parse_command["win"].format(dll=dll))
    elif os.name == "posix":
        # Linux
        os.system(parse_command["linux"].format(dll=dll))


if __name__ == '__main__':
    parse("api")
