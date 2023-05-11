import os
import re
import sys
from collections import defaultdict
from pathlib import Path
from typing import List, Tuple, Iterator, Dict

# Constant patterns
CPP_COMMENT_PATTERN = re.compile(r'/\*DOC_EXTRACT(.*?)\*/', re.DOTALL)
PYTHON_COMMENT_PATTERN = re.compile(r"'''(.*?)'''", re.DOTALL)
HEADER_PATTERN = re.compile(r'(##+.*?)\n', re.DOTALL)

def convert_inline_math(text: str) -> str:
    pattern = r"(?<!\$)(?<!\\)\$(?!\$)(?!`)(.*?)(?<!`)(?<!\\)\$(?!\$)"
    text = re.sub(pattern, r"$`\1`$", text)
    return text

def convert_math_star(text: str) -> str:
    pattern = r"((?<=\$`)(.*?)(?=`\$))|((?<=\$\$)(.*?)(?=\$\$))"
    text = re.sub(pattern, lambda m: m.group().replace("^*", "^\\ast"), text)
    return text

def highlight_keywords(text: str) -> str:
    text = convert_inline_math(text)
    text = convert_math_star(text)
    keyword_patterns = {
        'NOTE': (r'^NOTE:?\s*', '💡'),
        'WARNING': (r'^WARNING:?\s*', '⚠️'),
        'TODO': (r'^TODO:?\s*', '📝'),
        'IMPORTANT': (r'^IMPORTANT:?\s*', '❗'),
        'TIP': (r'^TIP:?\s*', '🌟'),
        'CHECKED': (r'^CHECKED:?\s*', '✅'),
    }

    for keyword, (pattern, emoji) in keyword_patterns.items():
        text = re.sub(pattern, f'**{emoji}**', text, flags=re.MULTILINE)

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

        # Extract the keyword at the top of the comment
        comment_parts = comment.split(maxsplit=1)
        keyword = comment_parts[0] if comment_parts else "DEFAULT"
        
        # Remove the keyword from the comment
        cleaned_comment = comment_parts[1] if len(comment_parts) > 1 else ""

        cleaned_comment = highlight_keywords(cleaned_comment)

        cleaned_comment = re.sub(r'!\[(.*?)\]\((.*?)\)', lambda m: f'![{m.group(1)}]({Path(input_file).parent / m.group(2)})', cleaned_comment)

        header_line = ""
        header_match = HEADER_PATTERN.search(cleaned_comment)
        if header_match:
            header_line = header_match.group(0)
            cleaned_comment = cleaned_comment.replace(header_line, '')
            headers_info.append((header_line.strip(), start_line))

        if header_line:
            keyword_comments[keyword].append(header_line)

        keyword_comments[keyword].append(cleaned_comment.strip() + '\n\n')
        keyword_comments[keyword].append(f'[{input_file}#L{start_line}]({input_file}#L{start_line})\n\n')
    
    return keyword_comments, headers_info

def generate_contents_table(headers_info: List[Tuple[str, int]], numbered: bool=False) -> str:
    contents_table = '# Contents\n\n'
    curr_section = 1
    curr_subsection = 0
    curr_subsubsection = 0
    for header, line_num in headers_info:
        if header.startswith("## "):
            prefix = f"{curr_section}. " if numbered else "- "
            contents_table += f"{prefix}[{header[3:]}](#{header[3:].replace(' ', '-')})\n"
            curr_section += 1
            curr_subsection = 0
        elif header.startswith("### "):
            curr_subsection += 1
            prefix = f"    {curr_section - 1}.{curr_subsection}. " if numbered else "    - "
            contents_table += f"{prefix}[{header[4:]}](#{header[4:].replace(' ', '-')})\n"
            curr_subsubsection = 0
        elif header.startswith("#### "):
            curr_subsubsection += 1
            prefix = f"        {curr_section - 1}.{curr_subsection}.{curr_subsubsection}. " if numbered else "        - "
            contents_table += f"{prefix}[{header[5:]}](#{header[5:].replace(' ', '-')})\n"
        contents_table += "\n"
    return contents_table + '\n\n'

def search_files(directory: str, extensions: Tuple[str, ...]) -> Iterator[str]:
    files = []
    for dirpath, _, filenames in os.walk(directory):
        for filename in filenames:
            if filename.endswith(extensions):
                files.append(os.path.join(dirpath, filename))
    return sorted(files)  # return sorted file paths

if __name__ == "__main__":
    if len(sys.argv) < 3:
        print("Usage: python extract_comments.py output_file.md search_directory")
        sys.exit(1)

    output_file = sys.argv[1]
    search_directory = sys.argv[2]
    file_extensions = (".cpp", ".hpp")

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

    contents_table = generate_contents_table(all_headers_info)

    with open(output_file, 'w') as md_file:
        md_file.write(contents_table)
        for keyword, comments in all_extracted_comments.items():
            # md_file.write(f"\n## {keyword}\n")
            md_file.write("\n".join(comments))
            md_file.write("\n---\n")  # Horizontal line after the entire group
        for comment in no_keyword_comments:
            md_file.write(comment)
            md_file.write("\n---\n")  # Horizontal line after each comment without a keyword
