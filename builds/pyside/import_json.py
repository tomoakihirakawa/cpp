import json
import asyncio
import websockets

async def server(websocket):
    data = await websocket.recv()

    data = json.loads(data)

    code = data.get('text')
    file = data.get('file')

    # Send the message back to the extension to be shown 
    # in the output window
    await websocket.send(code)

# host should always be localhost and
# port should match server<->client
start_server = websockets.serve(server, 'localhost', 54321)

asyncio.get_event_loop().run_until_complete(start_server)
asyncio.get_event_loop().run_forever()
