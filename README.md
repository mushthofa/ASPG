# ASPG
ASP-based attractor computation


PREREQUISITES:

1. Install Python 3 on your system
2. Download and compile the latest version of iclingo solver at: 

http://sourceforge.net/projects/potassco/files/iclingo/

and make sure that it is placed in the binary path of the system, i.e., you can call the program using:

(A source distribution is also included)

$ iclingo -v
iclingo 3.0.5 (clasp 1.3.10)
....

HOW TO USE:

1. Write a file to describe the network structure in ASP format. Each node is declared using the ASP fact:
	protein(X).
while each activation/inhibition between two nodes is declared using the ASP fact:
	activates(X, Y).
or
	inhibits(X, Y).

For example:

	% Proteins a,b and c where a activates b and c inhibits a
	protein(a).
	protein(b).
	protein(c).

	activates(a,b).
	inhibits(c,a).
	
2. Call the attractor computation program using the syntax:

	$ python3 attractor.py filename size [0|1] [a|s]

where filename is the file containing the network description, 0=activation rule r*, 1=activation rule r+, a=asynchronous update, s=synchronous update.
For example, if the network description above is stored in file 'network.txt', then to compute the attractors of the network when using activation rule r* and synchronous update, run the following command:

	$ python3 attractor.py network.txt 3 0 s
	
Some example network description files: 'toy-3nodes.txt', 'toy-4nodes.txt', 'yeast.txt' and 'th.txt' are provided.

