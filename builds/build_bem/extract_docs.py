import re
import sys
import markdown


def highlight_keywords(text):
    keyword_styles = {
        'NOTE': (r'(?i)(NOTE:?)', 'font-weight: bold; color: #1f77b4;'),
        'WARNING': (r'(?i)(WARNING:?)', 'font-weight: bold; color: #d62728;'),
        'TODO': (r'(?i)(TODO:?)', 'font-weight: bold; color: #ff7f0e;'),
        'IMPORTANT': (r'(?i)(IMPORTANT:?)', 'font-weight: bold; color: #9467bd;'),
        'TIP': (r'(?i)(TIP:?)', 'font-weight: bold; color: #2ca02c;'),
    }

    for keyword, (pattern, style) in keyword_styles.items():
        text = re.sub(
            pattern, f'<span style="{style}">{keyword}:</span>', text)

    return text


def extract_html_comments(input_file, output_file):
    markdown_comment_pattern = re.compile(r'/\*\*(.*?)\*/', re.DOTALL)

    with open(input_file, 'r') as cpp_file:
        content = cpp_file.read()

    markdown_comments = markdown_comment_pattern.finditer(content)

    with open(output_file, 'w') as html_file:
        # Add the CSS styles
        html_file.write('''
        <style>
        .note { font-weight: bold; color: #1f77b4; }
        .warning { font-weight: bold; color: #d62728; }
        .todo { font-weight: bold; color: #ff7f0e; }
        .important { font-weight: bold; color: #9467bd; }
        .tip { font-weight: bold; color: #2ca02c; }
        </style>
        ''')

        for match in markdown_comments:
            comment = match.group(1)
            start_line = content[:match.start()].count('\n') + 1

            # Remove leading asterisks and whitespace
            cleaned_comment = re.sub(
                r'^\s*\*', '', comment, flags=re.MULTILINE)

            cleaned_comment = highlight_keywords(cleaned_comment)

            # Convert the comment from Markdown to HTML
            html_comment = markdown.markdown(cleaned_comment)

            html_file.write(
                f'<p><a href="{input_file}#L{start_line}">[{input_file}#L{start_line}]</a>:</p>\n')
            html_file.write(html_comment + '\n')


if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("Usage: python extract_comments.py input_file.cpp output_file.html")
        sys.exit(1)

    input_file = sys.argv[1]
    output_file = sys.argv[2]
    extract_html_comments(input_file, output_file)
