from flask import Flask, render_template_string
from watchdog.observers import Observer
from watchdog.events import FileSystemEventHandler
import markdown2
import re
import threading
import logging

app = Flask(__name__)
markdown_html = ""
logging.basicConfig(level=logging.INFO)

class MyHandler(FileSystemEventHandler):
    def on_modified(self, event):
        global markdown_html
        if event.src_path.endswith('.cpp'):
            try:
                markdown_html = convert_file_to_html(event.src_path)
                logging.info(f"Updated HTML content from {event.src_path}")
            except Exception as e:
                logging.error(f"Error converting file to HTML: {e}")

def convert_file_to_html(file_path):
    try:
        with open(file_path, 'r') as file:
            content = file.read()
        comments = re.findall(r'/\*DOC_EXTRACT(.*?)\*/', content, re.DOTALL)
        markdown_text = '\n'.join(comments)
        return markdown2.markdown(markdown_text)
    except IOError as e:
        logging.error(f"Error reading file: {e}")
        return ""

@app.route('/')
def index():
    return render_template_string('<html><body>{{ content|safe }}</body></html>', content=markdown_html)

def start_flask_app():
    app.run(debug=True, use_reloader=False)

def start_watcher(path_to_watch):
    event_handler = MyHandler()
    observer = Observer()
    observer.schedule(event_handler, path=path_to_watch, recursive=False)
    observer.start()
    observer.join()

if __name__ == "__main__":
    file_to_watch = '/Users/tomoaki/Library/CloudStorage/Dropbox/code/cpp/builds/build_Fourier/example0_simple.cpp'  # Replace with your file path
    flask_thread = threading.Thread(target=start_flask_app)
    flask_thread.start()
    start_watcher(file_to_watch)