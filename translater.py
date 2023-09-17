import openai
import os
import sys

# Fetch the API key from environment variables
api_key = os.getenv("OPENAI_API_KEY")

# Raise an error if API key is not set
if api_key is None:
    raise ValueError("API Key must be set as an environment variable")

# Function to translate Japanese text to English
def translate_text(text):
    # Initialize the GPT-3 model
    model = "text-davinci-002"  # Replace with a valid model name

    # Translate Japanese text to English
    response = openai.Completion.create(
        engine=model,
        prompt=f"Translate the following Japanese document into English:\n'{text}'",
        max_tokens=100  # Adjust token limit as needed
    )

    # Fetch the translated text
    translated_text = response.choices[0].text.strip()

    return translated_text

if len(sys.argv) != 3:
    print("Usage: python script_name.py japanese.txt english.txt")
    sys.exit(1)

input_file = sys.argv[1]
output_file = sys.argv[2]

try:
    # Open file and read the content
    with open(input_file, 'r') as filehandle:
        japanese_text = filehandle.read()

    # Break the text into smaller chunks if it's too large
    max_chunk_size = 4000  # Make sure this is smaller than the model's maximum context size
    chunks = [japanese_text[i:i+max_chunk_size] for i in range(0, len(japanese_text), max_chunk_size)]

    translated_text_chunks = []
    for chunk in chunks:
        translated_text = translate_text(chunk)
        translated_text_chunks.append(translated_text)

    translated_text_final = ' '.join(translated_text_chunks)

    # Open file and write the content
    with open(output_file, 'w') as filehandle:
        filehandle.write(translated_text_final)

except FileNotFoundError:
    print(f"File {input_file} not found. Please make sure it exists.")
except Exception as e:
    print(f"An error occurred: {e}")
