#!/usr/bin/python
# -*- coding: utf-8 *-*
import re, os, shutil, sys

#===============================================================================
# Finds attractors in the network by looking for cycles
#===============================================================================
def attractors(filename, number_of_nodes, max_iterations):
  file = os.popen(
      "iclingo {filename} -c n={nodes} --imax={iterations} 0"
      .format(
          filename=filename, nodes=number_of_nodes, iterations=max_iterations
        )
    )
  
  # to a tuple (node, time, Bool) where Bool is True or False dependin
  # on whether the node is active or inhibited at time time
  a_list = re.findall(r'Answer: \d+\n(.*)\n', file.read())
  pat = re.compile(r'active\((\w+),(\d+)\)|inhibited\((\w+),(\d+)\)')
  a_list = [
    [
      y[:2] + (True,) if y[0] != '' else y[2:] + (False,)
      for y in re.findall(pat, x)
    ]
    for x in a_list
  ]

  a_sets = set()
  for x in a_list:
    a_set_t = tuple()
    for i in range(1, max_iterations + 1):
      t = tuple(sorted([(y[0], y[2]) for y in x if y[1] == str(i)]))
      if t:
        a_set_t = a_set_t + (t,)
    a_sets.add(tuple(sorted(a_set_t)))

  return a_sets

#===============================================================================
# Forbid attractors from showing themselves again
#===============================================================================
def add_visited_states(filename, visited):
  # Adds the visited states constraint to the ASP
  # stored in filename
  f = open(filename, 'a')
  for attractor in visited:
    for state in attractor:
      print(file=f), 
      print (
        "#cumulative t.\n"
        ":-"
        , file=f, end=" "),
      for gene, activated in state:
        if activated:
          print ("active({0},t),".format(gene), file=f, end=" "),
        else:
          print ("inhibited({0},t),".format(gene), file=f, end=" "),
      print ( "t > 0.", file=f)

      print (
        "#base.\n"
        ":- ", file=f, end=" "
        ),
      for gene, activated in state:
        if activated:
          print ( "active({0},0),".format(gene), file=f, end=" "),
        else:
          print ( "inhibited({0},0),".format(gene), file=f, end=" "),
      print ("t > 0.", file=f)

#===============================================================================
# Tells whether we can find a path of length max_time in the ASP program
#===============================================================================
def path_exists(filename, nodes, max_iterations):
  file = os.popen(
      "iclingo {filename} -c n={nodes} --iquery={time} | grep UNSATISFIABLE"
      .format(filename=filename, time=max_iterations, nodes=nodes)
    )

  return len(file.read()) == 0

  
#===============================================================================
# Write main ASP program
#===============================================================================
def write_main_sy(file_base, file_out, cycle, act):
    f = open(file_out, 'w')
    print(
    """
    #domain index_0(A).
    #domain index_0(I).


    #hide.
    #show active(X,T).
    #show inhibited(X,T).

    #const n=1.


    #base.
    index(1..n).
    index_0(0..n).
    
    %%%%%%%%%%%%%%%
    % Frame Axioms
    %%%%%%%%%%%%%%%

    % Activation and inhibition rules (A protein is activated if at
    % least one activator is on and no inhibitor is on)
    #cumulative t.
    """, file=f)
    
    if act == '+' :
        print("""
        active(Y,t)    :- protein(Y),
                  number_of_activation_links(Y,A,t-1),
                  number_of_inhibition_links(Y,I,t-1),
                  A>0, I=0, t>0.
                  
        inhibited(X,t) :- protein(X), not active(X,t).
        """, file=f)
    else:
        print("""
        active(Y,t)    :- protein(Y),
                  number_of_activation_links(Y,A,t-1),
                  number_of_inhibition_links(Y,I,t-1),
                  A>I, t>0.
                  
        inhibited(X,t) :- protein(X), not active(X,t).
        """, file=f)
    
    print("""
        % Counting links
        number_of_activation_links(Y,A,t) :- A {active(X,t) : activates(X,Y)} A, protein(Y).
        number_of_inhibition_links(Y,I,t) :- I {active(X,t) :  inhibits(X,Y)} I, protein(Y).

        %%%%%%%%%%%%%%%%%%%%%
        % Generate the cases (just initial conditions should be random)
        #base.

        number_of_activation_links(Y,A,0) :- A {active(X,0) :  activates(X,Y)} A, protein(Y).
        number_of_inhibition_links(Y,I,0) :- I {active(X,0) :  inhibits(X,Y)} I, protein(Y).

        active(X,0)    :- protein(X), not inhibited(X,0).
        inhibited(X,0) :- protein(X), not active(X,0).

        :- active(X,0), inhibited(X,0), protein(X).

        % Constraints
        #cumulative t.
        % A protein is either active or inhibited at any time point, i.e. the states are well defined
        :- active(X,t), inhibited(X,t), protein(X), t>=0.
        """, file=f)
        
    if cycle==True :
        print("""
        % Useful for finding all size 1 attractors
        #volatile t.
        :- not cycle(t).

        #cumulative t.
        part_of_cycle(Y, t) :- active(Y, t), active(Y, 0), protein(Y).
        part_of_cycle(Y, t) :- inhibited(Y, t), inhibited(Y, 0), protein(Y).

        cycle(t) :- n { part_of_cycle(Y, t): protein(Y) } n.   
        """, file=f)
        
    print("#base.", file=f)
    
    fbase = open(file_base, 'r')
    network = fbase.read()
    print(network, file=f)


  
#===============================================================================
# Write main ASP program
#===============================================================================
def write_main_as(file_base, file_out, cycle, act):
    f = open(file_out, 'w')
    print(
    """
    #domain index_0(A).
    #domain index_0(I).


    #hide.
    #show active(X,T).
    #show inhibited(X,T).

    #const n=1.


    #base.
    index(1..n).
    index_0(0..n).


    %%%%%%%%%%%%%%%
    % Frame Axioms
    %%%%%%%%%%%%%%%


    % Activation and inhibition rules (A protein is activated if at
    % least one activator is on and no inhibitor is on)
    #cumulative t.
    notmax(P) :- protein(P), protein(P1), P1!=P, P1>P.
    max(P) :- protein(P), not notmax(P).
    """, file=f)
    
    if act == '+' :
        print("""
        active(Y,t)    :- protein(Y), max(Y), 
                  number_of_activation_links(Y,A,t-1),
                  number_of_inhibition_links(Y,I,t-1),
                  A>0, I=0, t>0.
                  
        inhibited(X,t) :- protein(X), not active(X,t).
        """, file=f)
    else:
        print("""
        active(Y,t)    :- protein(Y), max(Y), 
                  number_of_activation_links(Y,A,t-1),
                  number_of_inhibition_links(Y,I,t-1),
                  A>I, t>0.
                  
        inhibited(X,t) :- protein(X), not active(X,t).
        """, file=f)
    
    print("""
        active(X, t) :- protein(X), not max(X), active(X, t-1), t>0.
        inhibited(X, t) :- protein(X), not max(X), inhibited(X, t-1), t>0.
        % Counting links
        number_of_activation_links(Y,A,t) :- A {active(X,t) : activates(X,Y)} A, protein(Y).
        number_of_inhibition_links(Y,I,t) :- I {active(X,t) :  inhibits(X,Y)} I, protein(Y).

        %%%%%%%%%%%%%%%%%%%%%
        % Generate the cases (just initial conditions should be random)
        #base.

        number_of_activation_links(Y,A,0) :- A {active(X,0) :  activates(X,Y)} A, protein(Y).
        number_of_inhibition_links(Y,I,0) :- I {active(X,0) :  inhibits(X,Y)} I, protein(Y).

        active(X,0)    :- protein(X), not inhibited(X,0).
        inhibited(X,0) :- protein(X), not active(X,0).

        :- active(X,0), inhibited(X,0), protein(X).

        % Constraints
        #cumulative t.
        % A protein is either active or inhibited at any time point, i.e. the states are well defined
        :- active(X,t), inhibited(X,t), protein(X), t>=0.
        %Prune states
        active(X, 1, X) :- protein(X), 
                        number_of_activation_links(X, A),
                        number_of_inhibition_links(X, I),
                        A>0, I=0.
                        
        inhibited(X, 1, X) :- protein(X), not active(X, 1, X).

        active(X, 1, C) :- protein(X),  protein(C), X!=C, active(X, 0).
        inhibited(X, 1, C) :- protein(X),  protein(C), X!=C, inhibited(X, 0).

        number_of_activation_links(X,A) :- A {active(Y,0) :  activates(Y,X)} A, protein(X).
        number_of_inhibition_links(X,I) :- I {active(Y,0) :  inhibits(Y,X)} I, protein(X).


        :- protein(C1), protein(C2), C1!=C2, active(X, 1, C1), inhibited(X, 1, C2).
        """, file=f)
        
    if cycle==True :
        print("""
        % Useful for finding all size 1 attractors
        #volatile t.
        :- not cycle(t).

        #cumulative t.
        part_of_cycle(Y, t) :- active(Y, t), active(Y, 0), protein(Y).
        part_of_cycle(Y, t) :- inhibited(Y, t), inhibited(Y, 0), protein(Y).

        cycle(t) :- n { part_of_cycle(Y, t): protein(Y) } n.   
        """, file=f)
    
    
    print("#base.", file=f)
    
    fbase = open(file_base, 'r')
    network = fbase.read()
    print(network, file=f)

    
    
    


if __name__ == '__main__':
  FILENAME = None
  FILENAME_CYCLE = None
  FILENAME_NO_CY = None
  TEMP_CYCLE     = ''
  TEMP_NO_CYCLE  = ''
  SIZE = 1

  args = sys.argv[1:]
  if len(args) >= 3:
    FILENAME = args[0]
    SIZE = int(args[1])
    if args[2] == '0':
        ACT = '*'
    else:
        ACT = '+'
    if len(args) >= 4 and args[3] == 'a':
        UPD = 'a'
    else:
        UPD = 's'
    
    FILENAME_CYCLE = FILENAME + "-cycle"
    FILENAME_NO_CY = FILENAME + "-no-cycle"
    
    if UPD == 's' :
        write_main_sy(FILENAME, FILENAME_CYCLE, True, ACT)
        write_main_sy(FILENAME, FILENAME_NO_CY, False, ACT)
    else:
        write_main_as(FILENAME, FILENAME_CYCLE, True, ACT)
        write_main_as(FILENAME, FILENAME_NO_CY, False, ACT)
    
    
    TEMP_CYCLE = FILENAME_CYCLE + ".tmp"
    TEMP_NO_CYCLE = FILENAME_NO_CY + ".tmp"

    shutil.copy(FILENAME_CYCLE, TEMP_CYCLE)
    shutil.copy(FILENAME_NO_CY, TEMP_NO_CYCLE)
  else:
    print("python3 attractor.py filename size [0|1] [a|s]", file=sys.stderr)
    sys.exit(1)

  max_i = 2**SIZE
  i = SIZE 

  attr = set()
  while path_exists(TEMP_NO_CYCLE, SIZE, i) and i<=max_i: 
    #print("i = {}, ".format(i), end=" ")
    attr_temp = attractors(TEMP_CYCLE, SIZE, i)
    attr_temp.discard(())
    attr = attr.union(attr_temp)
    add_visited_states(TEMP_CYCLE, attr_temp)
    add_visited_states(TEMP_NO_CYCLE, attr_temp)
    #print("found {} attractors".format(len(attr)))
    if len(attr_temp) == 0:
        i *= 2

  print("Found {} attractors".format(len(attr)))
  for i, x in enumerate(attr):
    print ("Attr #{} = {}".format(i + 1, x))
