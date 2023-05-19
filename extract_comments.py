import markdown
import os
import re
import sys
from collections import defaultdict
from pathlib import Path
from typing import List, Tuple, Iterator, Dict

# Constant patterns
CPP_COMMENT_PATTERN = re.compile(r'/\*DOC_EXTRACT(.*?)\*/', re.DOTALL)
PYTHON_COMMENT_PATTERN = re.compile(r"'''(.*?)'''", re.DOTALL)
HEADER_PATTERN = re.compile(r'(#+\s.*?)\n', re.DOTALL)

def search_labels(directory: str, extensions: Tuple[str, ...]) -> Dict[str, Tuple[str, int]]:
    label_pattern = re.compile(r'\\label\{(.*?)\}')
    labels = {}
    for dirpath, _, filenames in os.walk(directory):
        for filename in filenames:
            if filename.endswith(extensions):
                file_path = os.path.join(dirpath, filename)
                with open(file_path, 'r') as file:
                    content = file.read()
                for match in label_pattern.finditer(content):
                    label = match.group(1)
                    start_line = content[:match.start()].count('\n') + 1
                    labels[label] = (file_path, start_line)
    return labels


def convert_math_underscore(text: str) -> str:
    patterns = [
        r"(?<=\$\$)(.*?)(?=\$\$)",  # between $$
        r"(?<=\$\$\$)(.*?)(?=\$\$\$)",  # between $$$$
        r"(?<=`math\s)(.*?)(?=`)",  # between `math`
        r"(?<=\$)(.*?)(?=\$)"  # between single $
    ]
    for pattern in patterns:
        matches = list(re.finditer(pattern, text, re.DOTALL))
        for match in reversed(matches):
            start, end = match.span()
            # Replace underscores that do not have a space before them
            new_text = re.sub(r"(?<! )_", " _", match.group())
            text = text[:start] + new_text + text[end:]
    return text

def convert_inline_math(text: str) -> str:
    pattern = r"(?<!\$)(?<!\\)\$(?!\$)(?!`)(.*?)(?<!`)(?<!\\)\$(?!\$)"
    text = re.sub(pattern, r"$`\1`$", text)
    return text

def convert_math_star(text: str) -> str:
    pattern = r"((?<=\$`)(.*?)(?=`\$))|((?<=\$\$)(.*?)(?=\$\$))"
    text = re.sub(pattern, lambda m: m.group().replace("^*", "^\\ast"), text)
    return text

def highlight_keywords(text: str) -> str:
    text = convert_math_underscore(text)
    text = convert_inline_math(text)
    text = convert_math_star(text)
    keyword_patterns = {
        'NOTE': (r'^NOTE:?\s*', '💡'),
        'WARNING': (r'^WARNING:?\s*', '⚠️'),
        'TODO': (r'^TODO:?\s*', '📝'),
        'IMPORTANT': (r'^IMPORTANT:?\s*', '❗'),
        'TIP': (r'^TIP:?\s*', '🌟'),
        'CHECKED': (r'^CHECKED:?\s*', '✅'),
        'NOTIMPLEMENTED': (r'NOTIMPLEMENTED:?\s*', '☐'),
        'IMPLEMENTED': (r'IMPLEMENTED?\s*', '☑️'),
        '###': (r'^###:?\s*', '### ⚓️'),
        '##': (r'^## :?\s*', '## ⛵️'),
        '#': (r'^# :?\s*', '# 🐋'),
    }

    for keyword, (pattern, emoji) in keyword_patterns.items():
        text = re.sub(pattern, f'{emoji} ', text, flags=re.MULTILINE)


    def replace_label(match):
        label = match.group(1)
        text = match.group(2)
        file, line = labels.get(label, ('#', ''))
        return f'[{text}]({file}#L{line})' if file != '#' else f'[{text}](not found)'

    # text = re.sub(r'\\ref\{(.*?)\}\{(.*?)\}', replace_label, text)
    text = re.sub(r'\\ref\{(.*?)\}\{(.*?)\}', replace_label, text)

    return text

def extract_markdown_comments(input_file: str) -> Tuple[Dict[str, List[str]], List[Tuple[str, int]]]:
    with open(input_file, 'r') as file:
        content = file.read()

    markdown_comments = list(CPP_COMMENT_PATTERN.finditer(content)) + list(PYTHON_COMMENT_PATTERN.finditer(content))

    # Initialize dictionary to store comments based on keywords
    keyword_comments = defaultdict(list)
    headers_info = []

    for match in markdown_comments:
        comment = match.group(1)
        start_line = content[:match.start()].count('\n') + 1

        # Split comment into lines
        comment_lines = comment.split('\n')

        # Get the minimum number of leading spaces from all lines
        min_indent = min((len(re.match(r'^(\s*)', line).group(1)) for line in comment_lines if line.strip()), default=0)

        # Remove leading whitespace from all lines
        comment_lines = [line[min_indent:] for line in comment_lines]
        comment_lines = [line.lstrip() for line in comment_lines]


        # Try to extract "DOC_EXTRACT" keyword from the first line
        doc_extract_pattern = re.compile(r'^DOC_EXTRACT\s*(.*)')
        doc_extract_match = doc_extract_pattern.match(comment_lines[0])
        if doc_extract_match:
            keyword = doc_extract_match.group(1).strip() or 'DEFAULT'
            # Remove the first line from the comment
            comment_lines = comment_lines[1:]
        else:
            keyword = comment_lines[0].split()[0] if comment_lines and comment_lines[0].split() else 'DEFAULT'

        comment = '\n'.join(comment_lines)

        comment = comment.replace(keyword, '', 1) if keyword != 'DEFAULT' else comment

        cleaned_comment = highlight_keywords(comment)

        cleaned_comment = re.sub(r'!\[(.*?)\]\((.*?)\)', lambda m: f'![{m.group(1)}]({Path(input_file).parent / m.group(2)})', cleaned_comment)

        keyword_comments[keyword].append(cleaned_comment.strip() + '\n\n')
        # keyword_comments[keyword].append(f'[{input_file}#L{start_line}]({input_file}#L{start_line})\n\n')
        # keyword_comments[keyword].append('<p align="right">' + f'[{input_file}#L{start_line}]({input_file}#L{start_line})' + '</p>\n\n')
        keyword_comments[keyword].append('<p align="right"><a href="' + f'{input_file}#L{start_line}' + '">' + f'{input_file}#L{start_line}' + '</a></p>\n\n')


        # Extract header information for the contents table
        headers = re.findall(HEADER_PATTERN, cleaned_comment)
        for header in headers:
            headers_info.append((header.strip(), start_line))
    
    extracted_comments = []

    # Sort and group comments based on
    extracted_comments = []

    # Sort and group comments based on keyword
    for keyword in sorted(keyword_comments.keys()):
        extracted_comments.extend(keyword_comments[keyword])
        headers_info.append((keyword, len(keyword_comments[keyword])))

    return keyword_comments, headers_info


def generate_contents_table(headers_info: List[Tuple[str, int]], numbered: bool = False) -> str:    
    contents_table = '# Contents\n\n'
    curr_section = 1
    curr_subsection = 0
    curr_subsubsection = 0
    prefix = ''
    for header, line_num in headers_info:
        if header.startswith("# "):
            prefix = f"{curr_section}. " if numbered else "- "            
            added = f"{prefix}[{header[2:]}](#{header[2:].replace(' ', '-')})\n"
            contents_table += added
            curr_section += 1
            curr_subsection = 0
            print("added =",added)
        elif header.startswith("## "):
            curr_subsection += 1
            prefix = f"    {curr_section - 1}.{curr_subsection}. " if numbered else "    - "
            added = f"{prefix}[{header[3:]}](#{header[3:].replace(' ', '-')})\n"
            contents_table += added
            curr_subsubsection = 0
            print("added =",added)
        elif header.startswith("### "):
            curr_subsection += 1
            prefix = f"        {curr_section - 1}.{curr_subsection}. " if numbered else "        - "
            added = f"{prefix}[{header[4:]}](#{header[4:].replace(' ', '-')})\n"
            contents_table += added
            curr_subsubsection = 0
            print("added =",added)
        elif header.startswith("#### "):
            curr_subsubsection += 1            
            prefix = f"            {curr_section - 1}.{curr_subsection}.{curr_subsubsection}. " if numbered else "            - "
            added = f"{prefix}[{header[5:]}](#{header[5:].replace(' ', '-')})\n"
            contents_table += added
            print("added =",added)
        print("header =", header)
        print("prefix =", prefix)
        

    contents_table += "\n"

    return contents_table


def search_files(directory: str, extensions: Tuple[str, ...]) -> Iterator[str]:
    files = []
    for dirpath, _, filenames in os.walk(directory):
        for filename in filenames:
            if filename.endswith(extensions):
                files.append(os.path.join(dirpath, filename))
    return sorted(files)  # return sorted file paths


if __name__ == "__main__":
    if len(sys.argv) < 3:
        print("Usage: python extract_comments.py output_file search_directory")
        sys.exit(1)

    output_file = sys.argv[1]
    search_directory = sys.argv[2]
    file_extensions = (".cpp", ".hpp", "r.py")

    # Add the search_labels line here
    if len(sys.argv) >= 4:
        labels = search_labels(sys.argv[3], file_extensions)
    else:
        labels = search_labels(search_directory, file_extensions)

    all_extracted_comments = defaultdict(list)
    no_keyword_comments = []
    all_headers_info = []

    for input_file in sorted(search_files(search_directory, file_extensions)):  # Sort file paths alphabetically
        extracted_comments, headers_info = extract_markdown_comments(input_file)
        if extracted_comments:
            for keyword, comments in extracted_comments.items():
                if keyword == 'DEFAULT':
                    no_keyword_comments.extend(comments)
                else:
                    all_extracted_comments[keyword].extend(comments)
            all_headers_info.extend(headers_info)
    
    contents_table = generate_contents_table(all_headers_info) + "\n---\n"

    with open(output_file, 'w') as md_file:
        md_file.write(contents_table)
        for keyword, comments in all_extracted_comments.items():
            md_file.write("\n".join(comments))
            md_file.write("\n---\n")  # Horizontal line after the entire group
        for comment in no_keyword_comments:
            md_file.write(comment)
            md_file.write("\n---\n")  # Horizontal line after each comment without a keyword