import re
import sys


def highlight_keywords(text):
    keyword_patterns = {
        'note': r'(?i)(NOTE:?)',
        'warning': r'(?i)(WARNING:?)',
        'todo': r'(?i)(TODO:?)',
        'important': r'(?i)(IMPORTANT:?)',
        'tip': r'(?i)(TIP:?)',
    }

    for keyword, pattern in keyword_patterns.items():
        text = re.sub(
            pattern, f'<span class="{keyword}">{keyword.capitalize()}:</span>', text)

    return text


def extract_markdown_comments(input_file, output_file):
    markdown_comment_pattern = re.compile(r'/\*\*(.*?)\*/', re.DOTALL)

    with open(input_file, 'r') as cpp_file:
        content = cpp_file.read()

    markdown_comments = markdown_comment_pattern.finditer(content)

    with open(output_file, 'w') as md_file:
        # Add the link to the stylesheet at the beginning of the file
        md_file.write('<link rel="stylesheet" href="styles.css">\n\n')

        for match in markdown_comments:
            comment = match.group(1)
            start_line = content[:match.start()].count('\n') + 1

            # Remove leading asterisks and whitespace
            cleaned_comment = re.sub(
                r'^\s*\*', '', comment, flags=re.MULTILINE)

            cleaned_comment = highlight_keywords(cleaned_comment)

            md_file.write(
                f'[{input_file}#L{start_line}]({input_file}#L{start_line}):\n\n')
            md_file.write(cleaned_comment.strip() + '\n\n')


if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("Usage: python extract_comments.py input_file.cpp output_file.md")
        sys.exit(1)

    input_file = sys.argv[1]
    output_file = sys.argv[2]
    extract_markdown_comments(input_file, output_file)
