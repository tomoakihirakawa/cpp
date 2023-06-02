from typing import List, Tuple
import os
import re


def parse_headers(readme: str) -> List[Tuple[str, int]]:
    headers_info = []
    with open(readme, 'r') as f:
        lines = f.readlines()
        for i, line in enumerate(lines):
            if line.startswith("#"):
                headers_info.append((line.strip(), i + 1))
    return headers_info


def generate_contents_table(readme: str, headers_info: List[Tuple[str, int]], numbered: bool = False) -> str:
    contents_table = f"## [{os.path.basename(os.path.dirname(readme))}]({readme})\n\n"
    curr_section = 1
    curr_subsection = 0
    prefix = ''
    for header, line_num in headers_info:
        if header.startswith("# "):
            prefix = f"{curr_section}. " if numbered else "- "
            added = f"{prefix}[{header[2:]}]({readme}#{header[2:].replace(' ', '-')})\n"
            contents_table += added
            curr_section += 1
            curr_subsection = 0
        elif header.startswith("## "):
            curr_subsection += 1
            prefix = f"    {curr_section - 1}.{curr_subsection}. " if numbered else "    - "
            added = f"{prefix}[{header[3:]}]({readme}#{header[3:].replace(' ', '-')})\n"
            contents_table += added
    return contents_table


def generate_summary_readme(readme_files):
    with open("README.md", "w") as f:
        f.write("# Contents\n\n")

        for readme in readme_files:
            headers_info = parse_headers(readme)
            contents_table = generate_contents_table(readme, headers_info)
            f.write(contents_table)


if __name__ == "__main__":
    build_directories = os.listdir("./builds")
    directories = [os.path.join("./builds", build_dir)
                   for build_dir in build_directories]
    readme_files = []
    for directory in directories:
        readme_file = os.path.join(directory, "README.md")
        if os.path.exists(readme_file) and os.path.isfile(readme_file):
            readme_files.append(readme_file)
    generate_summary_readme(readme_files)
