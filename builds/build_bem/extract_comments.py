import re
import sys


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
    cpp_comment_pattern = re.compile(r'/\*\*(.*?)\*/', re.DOTALL)
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
        cleaned_comment = re.sub(r'^\s*\*', '', comment, flags=re.MULTILINE)

        cleaned_comment = highlight_keywords(cleaned_comment)

        extracted_comments += f'[{input_file}#L{start_line}]({input_file}#L{start_line}):\n\n'
        extracted_comments += cleaned_comment.strip() + '\n\n'

    return extracted_comments


if __name__ == "__main__":
    if len(sys.argv) < 3:
        print("Usage: python extract_comments.py output_file.md input_file1.cpp input_file2.cpp ...")
        sys.exit(1)

    output_file = sys.argv[1]
    input_files = sys.argv[2:]

    all_extracted_comments = ""

    for input_file in input_files:
        extracted_comments = extract_markdown_comments(input_file)
        if extracted_comments:
            file_link_name = f'## {input_file}\n\n'
            all_extracted_comments += file_link_name + extracted_comments

    with open(output_file, 'w') as md_file:
        md_file.write(all_extracted_comments)
