import plotly
import plotly.graph_objs as go
import json

lines = []
with open("cpu.csv") as f:
    lines = f.readlines()
    lines = lines[1:]
    lines = [l.split(',') for l in lines]

#print(lines)

first_time = int(lines[0][0])

times = [int(l[0]) - first_time for l in lines]
values = [float(l[7]) for l in lines]

plotly.offline.plot({
    "data": [go.Scatter(x=times, y=values)],
    "layout": go.Layout(title="CPU usage", 
      titlefont=dict(
        family='Helvetica',
        size=26
      ),
      yaxis=dict(
        range=[0, 100],
        title='% CPU Utilization'
      ),
      xaxis=dict(
        title='Time'
      ))
    }, auto_open=False)

with open("queue.json") as f:
    queue_data = json.load(f)

q_times = []
q_times.append(int(lines[0][0]))
q_values = []
q_values.append(0)

for val in queue_data:
    #print("queue data is {}".format(val))
    q_times.append(int(val[0]))
    q_values.append(int(val[1]))

q_times = [ t - first_time for t in q_times ]


plotly.offline.plot({
    "data": [go.Scatter(x=q_times, y=q_values)],
    "layout": go.Layout(title="Central Queue size", 
      titlefont=dict(
        family='Helvetica',
        size=26
      ),
      yaxis=dict(
        title='Queue size'
      ),
      xaxis=dict(
        title='Time'
      ))
    }, auto_open=False, filename='queue-plot.html')
