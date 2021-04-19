"""6.009 Lab 3 -- Circuit Solver."""

# NO IMPORTS ALLOWED!

# Uncomment below and comment/rename the solveLinear defined in this file to
# use the sample solveLinear function.
# Remember to comment it out before submitting!

#from solve_linear_sample import solveLinear
    
def substituteEquation(equation, substitutedVariable, substitutionEquation):
    """
        Note that implementing this function is optional. You might want to
        consider implementing it to break up your code into more managable
        chunks.
        
        Given:
            equation: An equation represented by a dictionary that maps the
                      variables or 1 to its coefficient in the equation.
                      E.g. {1: 2, 'x': 2, 'y': 3} represents 2 + 2x + 3y = 0.
            substitutedVariable: The variable to be substituted out of the
                                 equation.
            substitutionEquation: The substitution equation represented as a
                                  dictionary.
        Return:
            finalEquation: A dictionary representing the resulting equation
                           after the substitution is performed. 
    """
    sub_v = substitutionEquation[substitutedVariable]
    if substitutedVariable in equation:
        e_v = equation[substitutedVariable]
    else:
        e_v = 0
    
    diff = e_v/sub_v
#    new_eqn = {}

    #print(equation)
    #print(substitutionEquation)
    
    for i in equation:
        if i in substitutionEquation:
            #print(i)
#            new_eqn[i] = equation[i] - diff*substitutionEquation[i]
            equation[i] -= diff*substitutionEquation[i]
        else:
#            new_eqn[i] = equation[i]
            continue
    
    for i in substitutionEquation:
        if i not in equation:
#            new_eqn[i] = - diff*substitutionEquation[i]
            equation[i] = - diff*substitutionEquation[i]
        else:
            continue
    
#    key = list(new_eqn.keys())
#    l = []
    for i in list(equation.keys()):
        if equation[i] == 0: del equation[i]
    
#    for i in l:
#        del new_eqn[i]
    
    return equation
    

def solveLinear(variables, equations, result = None):
    """
        Given:
            variables: A set of strings or tuples representing the independent
                       variables. E.g. {'x', 'y', 'z'}
            equations: A list of linear equations where each equation is
                       represented as a dictionary. Each dictionary maps the
                       variables or 1 to its coefficient in the equation.
                       E.g. {1: 2, 'x': 2, 'y': 3} represents 2 + 2x + 3y = 0.
                       Note that all variables may not appear in all of the
                       equations. Moreover, you may assume that the equations
                       are independent.
        Return:
            result: A dictionary mapping each variable to the numerical value
            that solves the system of equations. Assume that there is exactly
            one solution. Some inaccuracy as typical from floating point
            computations will be acceptable.
    """
    # when only one equation left
    if result == None:
        result = {}
    
#    sub_eqn = {}
#    sub_v = 0
    
    if len(equations) == 1:
        left_v = next(iter(variables))
#        for i in variables:
#            left_v = i
        if 1 in equations[0]:
            const = equations[0][1]
            coff_v = equations[0][left_v]
            result[left_v] = -const/coff_v
        else:
            result[left_v] = 0
            
        return result
    
#    length = len(equations[0])
#    for i in equations:
#        if len(i) < length:
#            sub_eqn = i
#            length = len(i)
    
    # pick up the substitution eqn (the eqn with least number of vars)
#    sort_eqn = sorted(equations, key = len)
#    sub_eqn = sort_eqn[0]
    sub_eqn = min(equations, key = lambda eqn: len(eqn))
    
#    sub_v = max(sub_eqn, key = lambda v: sub_eqn[v])
    
#    coff = 0
#    for i in sub_eqn:
#        if i == 1:
#            break
#        elif abs(sub_eqn[i]) > coff:
#            sub_v = i
#            coff = abs(sub_eqn[i])
#        else:
#            continue
    
    # pick up the substituted variable (the variable with largest coefficient)
    sort_sub = sorted(sub_eqn, key = lambda i: abs(sub_eqn[i]), reverse=True)
    # use for loop will be faster
    
#    print(sort_sub)
    if sort_sub != []:
        if sort_sub[0] == 1:
            sub_v = sort_sub[1]
        else:
            sub_v = sort_sub[0]
#    else:
#        sub_v 
#    sub_v = max(abs([key for key in sub_eqn.keys() if type(key) != int]))
        
#    print('sub_var is :' + str(sub_v))
#    print('sub_eqn is:' + str(sub_eqn))
    
    # Build up the new systems with less equations and less variables (after applying substitution)
    new_eqns = []
    for i in equations:
        if i == sub_eqn:
            continue
        else:
            new_eqns.append(substituteEquation(i, sub_v, sub_eqn))
    
#    print(new_eqns)
    
#    new_var = set()
#    for i in variables:
#        if i == sub_v:
#            continue
#        else:
#            new_var.add(i)
    variables.remove(sub_v)
    
    # recursively substitution        
    new_sys = solveLinear(variables, new_eqns, result)
#    new_sys = solveLinear(new_var, new_eqns, result)
#    print('new_sys')
#    print(new_sys)
    
    
    # plug in the obtained value
    for i in new_sys:
        if i in sub_eqn:
            if 1 in sub_eqn:
                sub_eqn[1] += new_sys[i] * sub_eqn[i]
            else:
                sub_eqn[1] = new_sys[i] * sub_eqn[i]
                
            del sub_eqn[i]
    
    if 1 in sub_eqn:
        result[sub_v] = -sub_eqn[1]/sub_eqn[sub_v]
    else:
        result[sub_v] = 0
    
#    print('result')
#    print(result)
    
    return result
                
def solveCircuit(junctions, wires, resistances, voltages):
    """
        Given:
            junctions:  A set of junctions. Each junction is labeled by a string
                        or a tuple.
            wires:      A dictionary mapping a unique wire ID (a string or tuple)
                        to a tuple of two elements representing the starting and
                        ending junctions of the wire, respectively. The set of
                        wire IDs is disjoint from the set of junction labels.
                        Note that although electricity can flow in either
                        directions, each wire between a pair of junctions will
                        appear exactly once in the list. Moreover, the starting
                        and ending junctions are distinct.
            resitances: A dictionary mapping each unique wire ID to a numeric
                        value representing the magnitude of the resistance of
                        the wire in Ohms. This dictionary has the same set of
                        keys as the wires dictionary.
            voltages:   A dictionary mapping each unique wire ID to a numeric
                        value representing the voltage (EMF or potential
                        difference) of the battery connected along the wire in 
                        Volts. The positive terminal of the battery is next to
                        the ending junction (as defined in the wires dictionary)
                        if the voltage is positive whereas it is next to the 
                        starting junction otherwise. This dictionary also has
                        the same set of keys as the wires dictionary.
        Return:
            result: A dictionary mapping the label of each wire to the current
                    it carries. The labels must be the keys in the wires
                    dictionary and the current should be considered positive if
                    it is flowing from the starting junction to the ending
                    junction as specified in the wires dictionary.
    """

    current = {}
    eqns = []
    var = set()
    
    # eqns of KCL
    for i in junctions:
        var.add(i)    # str in junctions represents the voltage at that junction
        current[i] = []
        for w in wires:
            if i in wires[w]:
                current[i].append(w)
    
#    print(current)
    for i in current:
        eqn = {}
        for w in current[i]:
            var.add(w)     # key of wires represents the current through that wire
            if wires[w][0] == i:
                eqn[w] = -1
            if wires[w][1] == i:
                eqn[w] = 1
        eqns.append(eqn)
#    print(eqns)
    
#    print(len(eqns))
    
#    for i in range(len(new_eqns) - 1):
#        print('i: ' + str(i))
#        for j in range(i+1, len(new_eqns)):
#            print('j: ' + str(j))
#            if new_eqns[j].keys() != new_eqns[i].keys():
#                continue
#            else:
#                div = set()
#                div.add(ji/ii for ji,ii in zip(new_eqns[j], new_eqns[i]))
#                if len(div) == 1:
#                    dependent = True
#                    break
#    if dependent == True:
#        eqns.pop(0)

#    M = []
#    for eqn in eqns:
#        list_eqn = []
#        for w in wires:
#            if w in eqn:
#                list_eqn.append(eqn[w])
#            else:
#                list_eqn.append(0)
#        M.append(list_eqn)
#                
#    deter_M = det(M)
#    if deter_M == 0:
        
    # ensure the eqns are independent
    eqns.pop()
                    
    # ohm's law
    for w in wires:
        ohm = {}
        if resistances[w] != 0:
            ohm[wires[w][0]] = -1
            ohm[wires[w][1]] = 1
            ohm[w] = resistances[w]
        else:
            ohm[wires[w][0]] = -1
            ohm[wires[w][1]] = 1
        
        if voltages[w] != 0:
            ohm[1] = -voltages[w]
#        ohm[w] = -resistances[w]
        eqns.append(ohm)
    
    # set the grounded junction
    eqns.append({next(iter(junctions)): 1})
    
#    for eqn in eqns[:]:
#        for i in list(eqn.keys()):
#            if eqn[i] == 0:
#                del eqn[i]
    # solve linear equations
    
#    print(var)
#    print(eqns)
    
    soln_total = solveLinear(var, eqns)
    soln = {}
    for i in wires:
        if i in soln_total:
            soln[i] = soln_total[i]
    
    # sort obtained solution
#    soln_list = sorted(soln_dict.items(), key = lambda x: int(x[0]))
#    soln = {}
#    for i in soln_list:
#        soln[i[0]] = i[1]
        
    return soln

def findMaximumDeviationJunction(junctions, wires, resistances, voltages, currents):
    """
        Note that this part is completely optional and would not contribute to your grade.
        
        Given:
            junctions:  A set of junctions. Each junction is labeled by a
                        string or a tuple.
            wires:      A dictionary mapping a unique wire ID (a string or tuple)
                        to a tuple of two elements representing the starting and
                        ending junctions of the wire respectively. The set of
                        wire IDs is disjoint from the set of junction labels.
                        Note that although electricity can flow in either
                        direction, each wire between a pair of junctions will
                        appear exactly once in the list. Moreover, the starting
                        and ending junctions are distinct.
            resitances: A dictionary mapping each unique wire ID to a numeric
                        value representing the magnitude of the resistance of
                        the wire in Ohms. This dictionary has the same set of
                        keys as the wires dictionary.
            voltages:   A dictionary mapping each unique wire ID to a numeric
                        value representing the voltage (EMF or potential
                        difference) of the battery connected along the wire in
                        Volts. The positive terminal of the battery is next to
                        the ending junction (as defined in the wires dictionary)
                        if the voltage is positive whereas it is next to the
                        starting junction otherwise. This dictionary also has 
                        the same set of keys as the wires dictionary.
            currents:   A dictionary mapping each unique wire ID to a numeric
                        value representing the indicated current flowing along
                        the wire. The format is identical to that of the output 
                        of the previous function. Note that the values will not
                        necessarily be correct.
        Return:
            result: A junction with the maximum deviation from current
                    conservation. Note that any junction with maximal deviation
                    will be accepted.
    """
    kcl = {}
    for i in wires:
        for j in wires[i]:
            if j not in kcl:
                kcl[j] = [i]    # key is junction and value is wire
            else:
                kcl[j].append(i)
    
    eqn_junc = {}
    for i in kcl:   # for each junction
        eqn_junc[i] = {}
        for w in kcl[i]:     # for each wire
            if wires[w][0] == i:
                eqn_junc[i][w] = -1   # store coefficient of each wire w passing through each junction i 
            if wires[w][1] == i:
                eqn_junc[i][w] = 1
    
    max_junc = next(iter(junctions))
    max_dev = 0
    for j in eqn_junc:
        dev = 0
        for w in eqn_junc[j]:
            dev += eqn_junc[j][w]*currents[w]
        if abs(dev) > max_dev:
            max_dev = abs(dev)
            max_junc = j
    
    return max_junc
        
def findMaximumDeviationLoop(junctions, wires, resistances, voltages, currents):
    """
        Note that this part is completely optional and would not contribute to your grade.
        
        Given:
            junctions:  A set of junctions. Each junction is labeled by a string
                        or a tuple.
            wires:      A dictionary mapping a unique wire ID (a string or tuple)
                        to a tuple of two elements representing the starting and
                        ending junctions of the wire respectively. The set of
                        wire IDs is disjoint from the set of junction labels.
                        Note that although electricity can flow in either
                        directions, each wire between a pair of junctions will
                        appear exactly once in the list. Moreover, the starting
                        and ending junctions are distinct.
            resitances: A dictionary mapping each unique wire ID to a numeric
                        value representing the magnitude of the resistance of 
                        the wire in Ohms. This dictionary has the same set of
                        keys as the wires dictionary.
            voltages:   A dictionary mapping each unique wire ID to a numeric
                        value representing the voltage (EMF or potential
                        difference) of the battery connected along the wire in
                        Volts. The positive terminal of the battery is next to
                        the ending junction (as defined in the wires dictionary)
                        if the voltage is positive whereas it is next to the 
                        starting junction otherwise. This dictionary also has
                        the same set of keys as the wires dictionary.
            currents:   A dictionary mapping each unique wire ID to a numeric
                        value representing the indicated current flowing along
                        the wire. The format is identical to that of the output
                        of the previous function. Note that the values will not
                        necessarily be correct.
        Return:
            result: A list of wires IDs representing the edges along a loop with
                    maximal (additive) deviation from Kirchoff's loop law.
                    The wires should be in order along the cycle but the
                    starting node and the direction may be arbitrary.
    """
    neighbor = {i:set() for i in junctions}
    # build up graph of circuit (key: each junction; values: its neighbor junctions)
    for w in wires: 
        neighbor[wires[w][0]].add(wires[w][1])
        neighbor[wires[w][1]].add(wires[w][0])
    
    # backtrack function (for each junctions, setting it as starting point and search its neighbor in order; if next position)
    # is the starting point, save the loop
    
    # question: (1) repeating loops (abcda and bcdab)
    # question (2)  same starting point, different loops
    
    def back_junc(target, next_j, path = []):
        if next_j == target:
            return path
        
        neigh_j = neighbor[next_j]
        for n in neigh_j:
            if n not in path:
                path.append(n)
                back_junc(target, )
            
        

if __name__ == '__main__':
    # additional code here will be run only when lab.py is invoked directly
    # (not when imported from test.py), so this is a good place to put code
    # used for testing.
    pass
#    equation = {'w': 3, 'x': -1, 'y': 3, 'z': 10}
#    substitutedVariable = 'y'
#    substitutionEquation = {'w': 2, 'x': 1, 'y': 2, 1: 11}
#    result = substituteEquation(equation, substitutedVariable, substitutionEquation)
#    eqns = [{'x': 1, 'y': -1, 'z': 2, 1: 2},{'x': 1, 'y': 1, 'z': 3, 1: 1},{'x': 11, 'z': -13, 1: -20}]
#    variable = {'x','y','z'}
#    variable = {'0','1','2','3'}
#    eqns = [{'0': 1, '1': -1}, {'3': 1}, {'0': -1, '2': 1}, {'1': 1, '2': -1, '3': -1}]
#    result = solveLinear(variable, eqns)
    
#    junctions = {'A', 'B', 'C', 'D', 'E', 'F', 'G', 'H'}
#    wires = {('D', 'H'): ('D', 'H'), ('A', 'E'): ('A', 'E'),
#             ('C', 'D'): ('C', 'D'), ('E', 'F'): ('E', 'F'),
#             ('G', 'H'): ('G', 'H'), ('B', 'F'): ('B', 'F'),
#             ('C', 'G'): ('C', 'G'), ('F', 'G'): ('F', 'G'),
#             ('H', 'E'): ('H', 'E'), ('B', 'C'): ('B', 'C'),
#             ('D', 'A'): ('D', 'A'), ('A', 'B'): ('A', 'B'),
#             'battery': ('A', 'G')}
#    resistances = {('D', 'H'): 60, ('A', 'E'): 9, ('C', 'D'): 60, ('E', 'F'): 60,
#                   ('G', 'H'): 36, ('B', 'F'): 60, ('C', 'G'): 36, ('F', 'G'): 36,
#                   ('H', 'E'): 60, ('B', 'C'): 60, ('D', 'A'): 9, ('A', 'B'): 9,
#                   'battery': 0}
#    voltages = {('D', 'H'): 0, ('A', 'E'): 0, ('C', 'D'): 0, ('E', 'F'): 0,
#                ('G', 'H'): 0, ('B', 'F'): 0, ('C', 'G'): 0, ('F', 'G'): 0,
#                ('H', 'E'): 0, ('B', 'C'): 0, ('D', 'A'): 0, ('A', 'B'): 0,'battery': -300}
#    soln = solveCircuit(junctions, wires, resistances, voltages)