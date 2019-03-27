Github link:
https://github.com/mcbrider5002/level_4_project/tree/submission
(This branch is frozen to submission time.)

To my markers:

NB: The CLI isn't complete. Don't be surprised if some of the commands there don't work - having a user interface wasn't a huge priority for this project.

The only dependencies for this project are some Python packages (which should be provided with an Anaconda distribution) and some data.
I've added a hopefully exhaustive list to requirements.txt.
I've provided a subset of the data (the full data wouldn't fit under the submission limit) - if you really require the full data for some reason, my supervisor has it.

In order to run, navigate to the project folder. Then, main.py can be used to run the project.
To generate the visualisations of the scoring functions used: use 'python main.py experiment score'.
To plot the cyclosporin mass spectrum: use 'python main.py graph' (sadly the file cannot be provided as a parameter).
To run the tests: use 'python main.py tests'.

There are a few other things in there that you can play around with if you want.
(The CLI will prompt you with the options if you enter the name of a nonexistent command.)
However, I won't guarantee that it works, or will compute in reasonable time - the most important commands are up there.