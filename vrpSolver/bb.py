from gurobipy import * 
import numpy as np 
import copy 
import matplotlib.pyplot as plt 

class BBTreeNode:
    def __init__(self):
        self.lb = 0
        self.ub = float('inf')
        self.sol = {}
        self.intSol = {}
        self.branchVarList = []
        self.grbModel = None
        self.cnt = None
        self.isInteger = False



class BBTree:
    def __init__(self):
        self.root = BBTreeNodelNil()


RLP = Model('relaxed MIP')
x = {}
for i in range(14):
    x[i] = RLP.addVar(lb = 0 ,ub = GRB.INFINITY,vtype = GRB.CONTINUOUS, name = 'x_' + str(i))

RLP.setObjective(200 * x[0] + 3 * x[1] + 20 * x[2] + 50 * x[3] + 70 * x[4] + 20 * x[5] +
                 5 * x[6] + 10 * x[7] + 200 * x[8] + 150 * x[9] + 18 * x[10] + 8 * x[11] +
                 300 * x[12] + 185 * x[13], GRB.MAXIMIZE)

RLP.addConstr(6 * x[0] + 2 * x[1] + 3 * x[2] + x[6] + 4 * x[8] + 5 * x[11] <= 10, name = 'c_1')
RLP.addConstr(3 * x[1] + 5 * x[2] + 5 * x[4] + 8 * x[6] + 5 * x[8] + 8 * x[9]
              + 7 * x[11] + x[12] + 4 * x[13] <= 12, name = 'c_2')
RLP.addConstr(8 * x[4] + 1 * x[5] + 4 * x[9] + 2 * x[10] + 4 * x[12] + 5 * x[13] <= 14, name = 'c_3')
RLP.addConstr(8 * x[5] + 5 * x[7] + 7 * x[10] + x[12] + 3 * x[13] <= 14, name = 'c_4')
RLP.addConstr(10 * x[3] + 4 * x[5] + x[12] + 3 * x[13] <= 14, name = 'c_5')
RLP.addConstr(x[3] + x[4] <= 1, name = 'c_6')
RLP.addConstr(x[7] + x[10] <= 1, name = 'c_7')
RLP.addConstr(x[8] + x[13] <= 1, name = 'c_8')
RLP.addConstr(x[10] + x[1] <= 1, name = 'c_9')
RLP.addConstr(x[3] + x[2] <= 1, name = 'c_10')
RLP.addConstr(x[4] + x[2] <= 1, name = 'c_11')
RLP.addConstr(x[5] + x[2] <= 1, name = 'c_12')
RLP.addConstr(x[6] + x[2] <= 1, name = 'c_13')

RLP.optimize()

class Node:
    def __init__(self):
        self.local_LB = 0
        self.local_UB = np.inf
        self.x_sol = {}
        self.x_int_sol = {}
        self.branch_var_list = []
        self.model = None
        self.cnt = None 
        self.isInteger = False         

    def deepcopy_node(node):
        new_node = Node()
        new_node.local_LB = 0
        new_node.local_UB = np.inf
        new_node.x_sol = copy.deepcopy(node.x_sol)
        new_node.x_int_sol = copy.deepcopy(node.x_int_sol) 
        new_node.branch_var_list = []
        new_node.model = node.model.copy()
        new_node.cnt = node.cnt 
        new_node.isInteger = node.isInteger
        return new_node

def Branch_and_bound(RLP):
    # initialize the initial node
    RLP.optimize()
    global_UB = RLP.ObjVal
    global_LB = 0
    eps = 1e-3
    incumbent_node = None
    Gap = np.inf

    # create initial node # 初始化根节点（就是初始问题的松弛问题求解）
    Queue = []
    node = Node()
    node.local_LB = 0
    node.local_UB = global_UB
    node.model = RLP.copy()
    node.model.setParam("OutputFlag", 0)
    node.cnt = 0
    Queue.append(node)

    cnt = 0
    Global_UB_change = []
    Global_LB_change = []
    while (len(Queue) > 0 and global_UB - global_LB > eps):
        # select the current node
        current_node = Queue.pop()
        cnt += 1

        # solve the current model
        current_node.model.optimize()
        Solution_status = current_node.model.Status

        '''
        OPTIMAL = 2
        INFEASIBLE = 3
        UNBOUNDED = 5 
        '''

        # check whether the current solution is integer and execute prune step
        '''
        isInteger : mark whether the current solution is integer solution 
        isPruned  : mark whether the current solution is pruned 
        '''
        isInteger = True
        isPruned = False
        if (Solution_status == 2):
            for var in current_node.model.getVars():
                current_node.x_sol[var.varName] = var.x
                print(var.VarName, ' = ', var.x)

                # # round the solution to get an integer solution
                current_node.x_int_sol[var.varName] = (int)(var.x)  # round the solution to get an integer solution
                if (abs((int)(var.x) - var.x) >= eps):
                    isInteger = False
                    current_node.branch_var_list.append(var.VarName)  # to record the candidate branch variables

            # update the LB and UB
            if (isInteger == True):
                # For integer solution node, update the LB and UB
                current_node.isInteger = True
                current_node.local_LB = current_node.model.ObjVal
                current_node.local_UB = current_node.model.ObjVal
                if (current_node.local_LB > global_LB):
                    global_LB = current_node.local_LB
                    incumbent_node = Node.deepcopy_node(current_node)
            if (isInteger == False):
                # For integer solution node, update the LB and UB also
                current_node.isInteger = False
                current_node.local_UB = current_node.model.ObjVal
                current_node.local_LB = 0
                for var_name in current_node.x_int_sol.keys():
                    var = current_node.model.getVarByName(var_name)
                    current_node.local_LB += current_node.x_int_sol[var_name] * var.Obj
                if (current_node.local_LB > global_LB or (
                        current_node.local_LB == global_LB and current_node.isInteger == True)):
                    global_LB = current_node.local_LB
                    incumbent_node = Node.deepcopy_node(current_node)
                    incumbent_node.local_LB = current_node.local_LB
                    incumbent_node.local_UB = current_node.local_UB 

            # prune by optimility
            if (isInteger == True):
                isPruned = True

            # prune by bound
            if (isInteger == False and current_node.local_UB < global_LB):
                isPruned = True

            Gap = round(100 * (global_UB - global_LB) / global_LB, 2)
            print('\n ------------ \n', cnt, '\t Gap  = ', Gap, '  %')

        elif (Solution_status != 2):
            isInteger = False
            isPruned = True
            continue

        if (isPruned == False):
            # select the branch variable
            branch_var_name = current_node.branch_var_list[0]
            left_var_bound = (int)(current_node.x_sol[branch_var_name])
            right_var_bound = (int)(current_node.x_sol[branch_var_name]) + 1

            # create two child nodes
            left_node = Node.deepcopy_node(current_node)
            right_node = Node.deepcopy_node(current_node)

            # create left child node
            temp_var = left_node.model.getVarByName(branch_var_name)
            left_node.model.addConstr(temp_var <= left_var_bound, name='branch_left_' + str(cnt))
            left_node.model.setParam("OutputFlag", 0)
            left_node.model.update()
            cnt += 1
            left_node.cnt = cnt

            # create right child node
            temp_var = right_node.model.getVarByName(branch_var_name)
            right_node.model.addConstr(temp_var >= right_var_bound, name='branch_right_' + str(cnt))
            right_node.model.setParam("OutputFlag", 0)
            right_node.model.update()
            cnt += 1
            right_node.cnt = cnt

            Queue.append(left_node)
            Queue.append(right_node)
            
            # update the global UB, explore all the leaf nodes 
            temp_global_UB = 0 
            for node in Queue:
                node.model.optimize()
                if(node.model.status == 2):
                    if(node.model.ObjVal >= temp_global_UB):
                        temp_global_UB = node.model.ObjVal
            
            
            global_UB = temp_global_UB    
            Global_UB_change.append(global_UB)
            Global_LB_change.append(global_LB)                

    # all the nodes are explored, update the LB and UB 
    global_UB = global_LB 
    Gap = round(100 * (global_UB - global_LB) / global_LB, 2)
    Global_UB_change.append(global_UB)
    Global_LB_change.append(global_LB) 
    
    print('\n\n\n\n') 
    print('-----------------------------------------')
    print('         Branch and Bound terminates     ')
    print('         Optimal solution found          ')
    print('-----------------------------------------')
    print('\nFinal Gap  = ', Gap, '  %') 
    print('Optimal Solution:', incumbent_node.x_int_sol)
    print('Optimal Obj:', global_LB) 

    return incumbent_node, Gap, Global_UB_change, Global_LB_change
