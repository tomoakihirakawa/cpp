import os
import re
import sys
from pathlib import Path

def highlight_keywords(text):
    keyword_patterns = {
        'NOTE': (r'(?i)(NOTE:?)', '💡'),
        'WARNING': (r'(?i)(WARNING:?)', '⚠️'),
        'TODO': (r'(?i)(TODO:?)', '📝'),
        'IMPORTANT': (r'(?i)(IMPORTANT:?)', '❗'),
        'TIP': (r'(?i)(TIP:?)', '🌟'),
    }

    for keyword, (pattern, emoji) in keyword_patterns.items():
        text = re.sub(pattern, f'**{emoji} {keyword}:**', text)

    return text

def extract_markdown_comments(input_file):
    cpp_comment_pattern = re.compile(r'/\*\*EXPOSE(.*?)\*/', re.DOTALL)
    python_comment_pattern = re.compile(r"'''(.*?)'''", re.DOTALL)

    with open(input_file, 'r') as file:
        content = file.read()

    cpp_comments = cpp_comment_pattern.finditer(content)
    python_comments = python_comment_pattern.finditer(content)
    markdown_comments = list(cpp_comments) + list(python_comments)
    extracted_comments = ""

    for match in markdown_comments:
        comment = match.group(1)
        start_line = content[:match.start()].count('\n') + 1

        # Remove leading asterisks and whitespace
        cleaned_comment = comment
        # cleaned_comment = re.sub(r'^\s*\*', '', comment, flags=re.MULTILINE)
        # cleaned_comment = re.sub(r'\*([^*]+)\*', r'**\1**', cleaned_comment)

        cleaned_comment = highlight_keywords(cleaned_comment)

        # Wrap LaTeX equations with double dollar signs
        cleaned_comment = re.sub(r'(?<!\$)\$(?!\$)(.+?)\$(?!\$)', r'$$\1$$', cleaned_comment)

        header_line = ""
        header_match = re.search(r'## (.*?)\n', cleaned_comment)
        if header_match:
            header_line = header_match.group(0)
            cleaned_comment = cleaned_comment.replace(header_line, '')

        if header_line:
            extracted_comments += header_line
        extracted_comments += cleaned_comment.strip() + '\n\n'
        extracted_comments += f'[{input_file}#L{start_line}]({input_file}#L{start_line})\n\n'

    return extracted_comments



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
    # file_extensions = (".cpp", ".c", ".h", ".py")
    file_extensions = (".cpp", ".hpp")

    all_extracted_comments = ""

    for input_file in search_files(search_directory, file_extensions):
        extracted_comments = extract_markdown_comments(input_file)
        if extracted_comments:
            file_link_name = f'## {input_file}\n\n'
            all_extracted_comments += extracted_comments
            all_extracted_comments += '\n --- \n'
    
    with open(output_file, 'w') as md_file:
        md_file.write(all_extracted_comments)
