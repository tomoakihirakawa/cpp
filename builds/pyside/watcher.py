'''

```bash
python3 -m pip install watchdog
```

'''
import sys
import subprocess
from watchdog.observers import Observer
from watchdog.events import PatternMatchingEventHandler

class ReRunner(PatternMatchingEventHandler):
    def __init__(self, command):
        super().__init__(patterns=["*.py"], ignore_directories=True)
        self.command = command

    def on_any_event(self, event):
        print(f"Detected changes in {event.src_path}. Restarting...")
        self.restart()

    def restart(self):
        if hasattr(self, 'process'):
            self.process.terminate()
        self.process = subprocess.Popen(self.command, shell=True)

if __name__ == "__main__":
    path = sys.argv[1] if len(sys.argv) > 1 else "."
    command = sys.argv[2] if len(sys.argv) > 2 else "python your_pyside_script.py"

    event_handler = ReRunner(command)
    observer = Observer()
    observer.schedule(event_handler, path, recursive=True)
    observer.start()

    print("Watching for file changes. Press Ctrl+C to stop.")
    try:
        event_handler.restart()  # Start the application for the first time
        observer.join()
    except KeyboardInterrupt:
        observer.stop()
    observer.join()
