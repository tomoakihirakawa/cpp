from flask import Flask, request, render_template, redirect
import json
import os

app = Flask(__name__)

# Initialize with empty list
objects = []


@app.route("/", methods=["GET", "POST"])
def home():
    if request.method == "POST":
        # Add new object to the list
        name = request.form.get("name")
        obj_type = request.form.get("type")
        obj_file = request.form.get("objfile")

        if all([name, obj_type, obj_file]):
            new_object = {"name": name, "type": obj_type, "objfile": obj_file}
            objects.append(new_object)
            return redirect("/")
        else:
            return "Please fill all fields.", 400

    return render_template("home.html", objects=objects)


@app.route("/load", methods=["POST"])
def load():
    settings_file = request.files.get('file')

    if settings_file:
        settings = json.load(settings_file)

        for obj in settings.get('input_files', []):
            if 'name' in obj and 'type' in obj and 'objfile' in obj:
                objects.append(obj)
            else:
                return f'Object missing "name", "type", or "objfile"', 400

        return redirect("/")
    else:
        return "No file uploaded.", 400


if __name__ == "__main__":
    app.run(debug=True)
