# # python3 extract_comments.py README.md ./

# import os
# import re
# import sys
# from pathlib import Path

# def highlight_keywords(text):
#     keyword_patterns = {
#         'NOTE': (r'(?i)(NOTE:?)', '💡'),
#         'WARNING': (r'(?i)(WARNING:?)', '⚠️'),
#         'TODO': (r'(?i)(TODO:?)', '📝'),
#         'IMPORTANT': (r'(?i)(IMPORTANT:?)', '❗'),
#         'TIP': (r'(?i)(TIP:?)', '🌟'),
#     }

#     for keyword, (pattern, emoji) in keyword_patterns.items():
#         text = re.sub(pattern, f'**{emoji} {keyword}:**', text)

#     return text

# def insert_space(match):
#     return ' ' + match.group(0) if match.group(1) != '\n' else match.group(0)

# def extract_markdown_comments(input_file):
#     cpp_comment_pattern = re.compile(r'/\*DOC_EXTRACT(.*?)\*/', re.DOTALL)
#     python_comment_pattern = re.compile(r"'''(.*?)'''", re.DOTALL)

#     with open(input_file, 'r') as file:
#         content = file.read()

#     cpp_comments = cpp_comment_pattern.finditer(content)
#     python_comments = python_comment_pattern.finditer(content)
#     markdown_comments = list(cpp_comments) + list(python_comments)
#     extracted_comments = ""

#     for match in markdown_comments:
#         comment = match.group(1)
#         start_line = content[:match.start()].count('\n') + 1

#         # Remove leading asterisks and whitespace
#         cleaned_comment = comment

#         cleaned_comment = highlight_keywords(cleaned_comment)

#         # Replace the relative image path with the correct path
#         cleaned_comment = re.sub(r'!\[(.*?)\]\((.*?)\)', lambda m: f'![{m.group(1)}]({os.path.join(os.path.dirname(input_file), m.group(2))})', cleaned_comment)

#         header_line = ""

#         header_match = re.search(r'## (.*?)\n', cleaned_comment)
#         if header_match:
#             header_line = header_match.group(0)
#             cleaned_comment = cleaned_comment.replace(header_line, '')

#         if header_line:
#             extracted_comments += header_line
#         extracted_comments += cleaned_comment.strip() + '\n\n'
#         extracted_comments += f'[{input_file}#L{start_line}]({input_file}#L{start_line})\n\n'

#     return extracted_comments


# def search_files(directory, extensions):
#     for dirpath, _, filenames in os.walk(directory):
#         for filename in filenames:
#             if filename.endswith(extensions):
#                 yield os.path.join(dirpath, filename)


# if __name__ == "__main__":
#     if len(sys.argv) < 3:
#         print("Usage: python extract_comments.py output_file.md search_directory")
#         sys.exit(1)

#     output_file = sys.argv[1]
#     search_directory = sys.argv[2]
#     # file_extensions = (".cpp", ".c", ".h", ".py")
#     file_extensions = (".cpp", ".hpp")

#     all_extracted_comments = ""

#     for input_file in search_files(search_directory, file_extensions):
#         extracted_comments = extract_markdown_comments(input_file)
#         if extracted_comments:
#             file_link_name = f'## {input_file}\n\n'
#             all_extracted_comments += extracted_comments
#             all_extracted_comments += '\n --- \n'
    
#     with open(output_file, 'w') as md_file:
#         md_file.write(all_extracted_comments)





# python3 extract_comments.py README.md ./

import os
import re
import sys
from collections import defaultdict
from pathlib import Path

def convert_inline_math(text):
    # This pattern will find all instances of $math expression$ but not $`math expression`$, $$math equation$$, and \$somthing\$
    pattern = r"(?<!\$)(?<!\\)\$(?!\$)(?!`)(.*?)(?<!`)(?<!\\)\$(?!\$)"
    # Replace matched patterns with new format
    text = re.sub(pattern, r"$`\1`$", text)
    return text

def highlight_keywords(text):
    text = convert_inline_math(text)  # Add this line to call the function before other replacements
    keyword_patterns = {
        'NOTE': (r'(?i)(NOTE:?)', '💡'),
        'WARNING': (r'(?i)(WARNING:?)', '⚠️'),
        'TODO': (r'(?i)(TODO:?)', '📝'),
        'IMPORTANT': (r'(?i)(IMPORTANT:?)', '❗'),
        'TIP': (r'(?i)(TIP:?)', '🌟'),
        'CHECKED': (r'(?i)(CHECKED:?)', '✅'),
    }

    for keyword, (pattern, emoji) in keyword_patterns.items():
        text = re.sub(pattern, f'**{emoji} {keyword}:**', text)

    return text

def insert_space(match):
    return ' ' + match.group(0) if match.group(1) != '\n' else match.group(0)

def generate_contents_table(headers_info, numbered=False):
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



def extract_markdown_comments(input_file):
    cpp_comment_pattern = re.compile(r'/\*DOC_EXTRACT(.*?)\*/', re.DOTALL)
    python_comment_pattern = re.compile(r"'''(.*?)'''", re.DOTALL)

    with open(input_file, 'r') as file:
        content = file.read()

    cpp_comments = cpp_comment_pattern.finditer(content)
    python_comments = python_comment_pattern.finditer(content)
    markdown_comments = list(cpp_comments) + list(python_comments)
    extracted_comments = ""
    headers_info = []

    for match in markdown_comments:
        comment = match.group(1)
        start_line = content[:match.start()].count('\n') + 1

        # Remove leading asterisks and whitespace
        cleaned_comment = comment

        cleaned_comment = highlight_keywords(cleaned_comment)

        # Replace the relative image path with the correct path
        cleaned_comment = re.sub(r'!\[(.*?)\]\((.*?)\)', lambda m: f'![{m.group(1)}]({os.path.join(os.path.dirname(input_file), m.group(2))})', cleaned_comment)

        header_line = ""
        header_match = re.search(r'(##+.*?)\n', cleaned_comment)
        if header_match:
            header_line = header_match.group(0)
            cleaned_comment = cleaned_comment.replace(header_line, '')
            headers_info.append((header_line.strip(), start_line))

        if header_line:
            extracted_comments += header_line
        extracted_comments += cleaned_comment.strip() + '\n\n'
        extracted_comments += f'[{input_file}#L{start_line}]({input_file}#L{start_line})\n\n'

    return extracted_comments, headers_info

def search_files(directory, extensions):
    for dirpath, _, filenames in os.walk(directory):
        for filename in filenames:
            if filename.endswith(extensions):
                yield os.path.join(dirpath, filename)


if __name__ == "__main__":
    if len(sys.argv) < 3:
        print("Usage: python extract_comments.py output_file.md search_directory")
        sys.exit(1)

    output_file = sys.argv[1]
    search_directory = sys.argv[2]
    file_extensions = (".cpp", ".hpp")

    all_extracted_comments = ""
    all_headers_info = []

    for input_file in search_files(search_directory, file_extensions):
        extracted_comments, headers_info = extract_markdown_comments(input_file)
        if extracted_comments:
            all_extracted_comments += extracted_comments
            all_extracted_comments += '\n --- \n'
            all_headers_info.extend(headers_info)

    contents_table = generate_contents_table(all_headers_info)

    with open(output_file, 'w') as md_file:
        md_file.write(contents_table)
        md_file.write(all_extracted_comments)
