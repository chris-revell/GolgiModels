
# Golgi Models

To use with Docker, from the cloned repository directory,
```
docker build -t golgimodels:latest .
```
This may take about 15 minutes.

To run,
```
docker run -it --rm -p 9384:9384 golgimodels:latest
```

It may take about 5 minutes before you get the `julia>` prompt,
indicating that the server is running.

Point your browser to http://localhost:9384/browser-display.  You
should see some plots and buttons.
