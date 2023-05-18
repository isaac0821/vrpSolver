from .common import *
from .error import *

class BSTreeNode(object):
    def __init__(self, key, value, parent=None, left=None, right=None):
        self.key = key
        self.value = value
        self.parent = parent
        self.left = left
        self.right = right

    def print(self):
        print("Key: ", self.key, 
            "Value: ", self.value, 
            "Parent: ", self.parent.key if self.parent != None else 'None',
            "Left: ", self.left.key if self.left != None else 'None',
            "Right: ", self.right.key if self.right != None else 'None')
        return

class BSTree(object):
    # Binary Search Tree ======================================================
    # NOTE: Father/ancestor of a series of tree structures

    # Initialize with an empty tree
    def __init__(self):
        self.nil = None
        self.root = self.nil

    # Query using key
    def query(self, key):
        def search(n, key):
            if (n == self.nil):
                return self.nil
            if (n.key == key):
                return n
            if (key < n.key):
                return search(n.left, key)
            else:
                return search(n.right, key)
        return search(self.root, key)

    def print(self):
        t = self.traverse()
        for n in t:
            n.print()
        return

    def traverse(self, mode='Left'):
        traverse = []
        def leftTraverse(n):
            if (n != self.nil):
                leftTraverse(n.left)
                traverse.append(n)
                leftTraverse(n.right)
            return
        def midTraverse(n):
            if (n != self.nil):
                traverse.append(n)
                midTraverse(n.left)                
                midTraverse(n.right)
            return
        def rightTraverse(n):
            if (n != self.nil):
                rightTraverse(n.right)
                traverse.append(n)
                rightTraverse(n.left)
            return
        if (mode == 'Left'):
            leftTraverse(self.root)
        elif (mode == 'Mid'):
            midTraverse(self.root)
        elif (mode == 'Right'):
            rightTraverse(self.root)
        return traverse

    # Insert a BSTreeNode to BST
    def insert(self, n):
        x = self.root
        y = self.nil
        while (x != self.nil):
            y = x
            if (n.key < x.key):
                x = x.left
            elif(n.key > x.key):
                x = x.right
            else:
                raise KeyExistError("Key %s exists in BST." % n.key)

        n.parent = y
        if (y == self.nil):
            self.root = n
        elif (n.key < y.key):
            y.left = n
        else:
            y.right = n
        return

    # Delete a BSTreeNode from BST
    def deleteByKey(self, key):
        n = self.query(key)
        if (n != self.nil):
            return self.delete(n)
        else:
            raise KeyNotExistError("Cannot find key %s in BST" % key)
    def delete(self, n):
        # Replace node u with node v in this previous location
        def replace(u, v):
            if (u.parent == self.nil):
                self.root = u
            elif (u == u.parent.left):
                u.parent.left = v
            else:
                u.parent.right = v
            if (v != self.nil):
                v.parent = u.parent

        if (n.left == self.nil):
            replace(n, n.right)
        elif (n.right == self.nil):
            replace(n, n.left)
        else:
            y = self.min(n.right)
            if (y != n.right):
                replace(y, y.right)
                y.right = n.right
                y.right.parent = y
            replace(n, y)
            y.left = n.left
            y.left.parent = y
        return

    # Predecessor
    def prevByKey(self, key):
        n = self.query(key)
        if (n != self.nil):
            return self.prev(n)
        else:
            raise KeyNotExistError("Cannot find key %s in BST" % key)
    def prev(self, n):
        if (n.left != self.nil):
            return self.max(n.left)
        else:
            y = n.parent
            while (y != self.nil and n == y.left):
                n = y
                y = y.parent
            return y

    # Successor
    def nextByKey(self, key):
        n = self.query(key)
        if (n != self.nil):
            return self.next(n)
        else:
            raise KeyNotExistError("Cannot find key %s in BST" % key)
    def next(self, n):
        if (n.right != self.nil):
            return self.min(n.right)
        else:
            y = n.parent
            while (y != self.nil and n == y.right):
                n = y
                y = y.parent
            return y

    # Get minimum
    def min(self, n=None):
        if (n == None):
            n = self.root
        if (n == self.nil):
            n = self.root
        while (n.left != self.nil):
            n = n.left
        return n

    # Get maximum
    def max(self, n=None):
        if (n == None):
            n = self.root
        if (n == self.nil):
            n = self.root
        while (n.right != self.nil):
            n = n.right
        return n

class RedBlackTreeNode(BSTreeNode):
    def __init__(self, key, value, parent=None, left=None, right=None, color=None):
        super(RedBlackTreeNode, self).__init__(key, value, parent, left, right)
        self.color = color

    def print(self):
        print("Key: ", self.key, 
            "Value: ", self.value, 
            "Color: ", self.color,
            "Parent: ", self.parent.key if self.parent != None else 'None',
            "Left: ", self.left.key if self.left != None else 'None',
            "Right: ", self.right.key if self.right != None else 'None')
        return

class RedBlackTree(BSTree):
    def __init__(self):
        self.nil = RedBlackTreeNode(key=None, value=None, color='B')
        self.root = self.nil        

    # Left rotation
    def leftRotation(self, n):
        y = n.right
        n.right = y.left
        if (y.left != self.nil):
            y.left.parent = n
        y.parent = n.parent
        if (n.parent == self.nil):
            self.root = y
        elif (n == n.parent.left):
            n.parent.left = y
        else:
            n.parent.right = y
        y.left = n
        n.parent = y
        return

    # Right rotation
    def rightRotation(self, n):
        y = n.left
        n.left = y.right
        if (y.right != self.nil):
            y.right.parent = n
        y.parent = n.parent
        if (n.parent == self.nil):
            self.root = y
        elif (n == n.parent.right):
            n.parent.right = y
        else:
            n.parent.left = y
        y.right = n
        n.parent = y
        return

    # Insertion
    def insert(self, n):
        x = self.root
        y = self.nil
        while (x != self.nil):
            y = x
            if (n.key < x.key):
                x = x.left
            elif(n.key > x.key):
                x = x.right
            else:
                raise KeyExistError("Key %s exists in BST." % n.key)
        n.parent = y
        if (y == self.nil):
            self.root = n
        elif (n.key < y.key):
            y.left = n
        else:
            y.right = n
        n.left = self.nil
        n.right = self.nil
        n.color = 'R'
        while (n.parent.color == 'R'):
            if (n.parent == n.parent.parent.left):
                y = n.parent.parent.right
                if (y.color == 'R'):
                    n.parent.color ='B'
                    y.color = 'B'
                    n.parent.parent.color = 'R'
                    n = n.parent.parent
                else:
                    if (n == n.parent.right):
                        n = n.parent
                        self.leftRotation(n)
                    n.parent.color = 'B'
                    n.parent.parent.color = 'R'
                    self.rightRotation(n.parent.parent)
            else:
                y = n.parent.parent.left
                if (y.color == 'R'):
                    n.parent.color ='B'
                    y.color = 'B'
                    n.parent.parent.color = 'R'
                    n = n.parent.parent
                else:
                    if (n == n.parent.left):
                        n = n.parent
                        self.rightRotation(n)
                    n.parent.color = 'B'
                    n.parent.parent.color = 'R'
                    self.leftRotation(n.parent.parent)
        self.root.color = 'B'
        return

    # Delete
    def delete(self, n):
        def replace(u, v):
            if (u.parent == self.nil):
                self.root = u
            elif (u == u.parent.left):
                u.parent.left = v
            else:
                u.parent.right = v
            if (v != self.nil):
                v.parent = u.parent
        y = n
        yOriColor = y.color
        if (n.left == self.nil):
            x = n.right
            replace(n, n.right)
        elif (n.right == self.nil):
            x = n.left
            replace(n, n.left)
        else:
            y = self.min(n.right)
            x = y.right
            if (y != n.right):
                replace(y, y.right)
                y.right = n.right
                y.right.parent = y
            else:
                x.parent = y
            replace(n, y)
            y.left = n.left
            y.left.parent = y
            y.color = n.color
        if (yOriColor == 'B'):
            while (x != self.root and x.color == 'B'):
                if (x == x.parent.left):
                    w = x.parent.right
                    if (w.color == 'R'):
                        w.color == 'B'
                        x.parent.color = 'R'
                        leftRotation(x.parent)
                        w = x.parent.right
                    if (w.left.color == 'B' and w.right.color == 'B'):
                        w.color = 'R'
                        x = x.parent
                    else:
                        if (w.right.color == 'B'):
                            w.left.color = 'B'
                            w.color = 'R'
                            rightRotation(w)
                            w = x.parent.right
                        w.color = x.parent.color
                        x.parent.color = 'B'
                        w.right.color = 'B'
                        leftRotation(x.parent)
                        x = self.root
                else:
                    w = x.parent.left
                    if (w.color == 'R'):
                        w.color == 'B'
                        x.parent.color = 'R'
                        rightRotation(x.parent)
                        w = x.parent.left
                    if (w.right.color == 'B' and w.left.color == 'B'):
                        w.color = 'R'
                        x = x.parent
                    else:
                        if (w.left.color == 'B'):
                            w.right.color = 'B'
                            w.color = 'R'
                            leftRotation(w)
                            w = x.parent.left
                        w.color = x.parent.color
                        x.parent.color = 'B'
                        w.left.color = 'B'
                        rightRotation(x.parent)
                        x = self.root
            x.color = 'B'
        return

class IntervalTreeNode(RedBlackTreeNode):
    def __init__(self, lower, upper, value, parent=None, left=None, right=None, color=None, childUpper=None):
        self.key = lower
        self.value = value
        self.lower = lower # float('-inf') represents infinite
        self.upper = upper # float('inf') represents infinite
        self.parent = parent
        self.left = left
        self.right = right
        self.color = color  # RedBlackTree

class NonOverlapIntervalTree(RedBlackTree):
    def __init__(self):
        self.nil = IntervalTreeNode(None, None, None, color='B')
        self.root = self.nil

    @property    
    def isEmpty(self):
        if (self.root == self.nil):
            return True
        else:
            return False

    def query(self, t):
        def search(n, t):
            if (n == self.nil):
                return self.nil
            if (n.lower <= t and t <= n.upper):
                return n
            if (t < n.lower):
                return search(n.left, t)
            elif (t > n.upper):
                return search(n.right, t)
        return search(self.root, t)

    def earlisetAvail(self, t):
        visited = []
        def search(n, t):
            if (n == self.nil):
                return self.nil
            else:
                visited.append(n)
            if (n.lower <= t and t <= n.upper):
                return n
            if (t < n.lower):
                return search(n.left, t)
            elif (t > n.upper):
                return search(n.right, t)
        s = search(self.root, t)
        if (s != self.nil):
            return t, s
        elif (t > visited[-1].upper):
            i = self.next(visited[-1])
            return i.lower, i
        elif (t < visited[-1].lower):
            return visited[-1].lower, visited[-1]

class VisitNode(object):
    def __init__(self, key, serviceTime=0):
        self.key = key
        self.serviceTime = serviceTime
        self.arr = 0
        self.dep = serviceTime
        self.prev = None
        self.succ = None

class RoundTripSym(object):
    # Linked list =============================================================
    # NOTE: Symmetric round trip for TSP/VRP, starts from and returns to the depot
    def __init__(self, tau, depotID, serviceTime = 0):
        # Initialize the round trip with a depot
        depot = VisitNode(depotID)
        depot.arr = 0
        depot.dep = depot.arr + depot.serviceTime
        self.head = depot
        self.tail = depot
        # Given (key, key) returns a distance
        self.tau = {}

    @property
    def length(self):
        return self.tail.dep

    @property
    def stops(self):
        return len([i for i in self.traverse()])

    def prepend(self, n):
        n.succ = self.head
        n.prev = None
        if (self.head != None):
            self.head.prev = n
        self.head = n
        return

    def append(self, n):
        if (self.head == None):
            self.head = n
            self.tail = n
        else:
            self.tail.succ = n
            n.prev = self.tail
            self.tail = n
        return

    def insertAfterKey(self, n, key):
        x = self.query(key)
        self.insert(n, x)
        return
    def insert(self, n, x):
        # Insert node n after x
        n.succ = x.succ
        n.prev = x
        if (x.succ != None):
            x.succ.prev = n
        x.succ = n
        return

    def deleteByKey(self, key):
        n = self.query(key)
        self.delete(n)
        return
    def delete(self, n):
        if (n.prev != None):
            n.prev.succ = n.succ
        else:
            self.head = n.succ
        if (n.succ != None):
            n.succ.prev = n.prev
        return

    def query(self, key):
        n = self.head
        while (n != None and n.key != key):
            n = n.succ
        return n

    def traverse(self):
        traverse = []
        n = self.head
        while (n.succ != None):
            traverse.append(n)
            n = n.succ
        traverse.append(n)
        return traverse



class JobNode(object):
    def __init__(self, key, value=None, ts=None, te=None, prev=None, succ=None):
        self.key = key
        self.ts = ts
        self.te = te
        self.prev = prev
        self.succ = succ

class JobLinkedList(object):
    def __init__(self, tws):
        self.tws = tws           # Time windows as a NonOverlapIntervalTree
        self.head = None
        self.tail = None         # tail要单独处理, 只有一个元素的时候tail是head的复制
        
    @property
    def makespan(self):
        if (self.tail != None):
            return self.tail.te
        else:
            return 0

    @property
    def length(self):
        l = len([i for i in self.traverse()])
        return l

    def query(self, key):
        n = self.head
        while (n != None and n.key != key):
            n = n.succ
        if (n == None):
            raise KeyNotExistError("%s does not exist" % key)  
        return n

    def clone(self):
        c = JobLinkedList(self.tws)
        tr = self.traverse()
        tc = []
        for i in range(len(tr)):
            tc.append(JobNode(
                key = tr[i].key,
                value = tr[i].key,
                ts = tr[i].ts,
                te = tr[i].te,
                prev = tc[i - 1] if i >= 1 else None))
        for i in range(len(tr)):
            tc[i].succ = tc[i + 1] if (i < len(tr) - 1) else None
        if (len(tc) > 0):
            c.head = tc[0]
            c.tail = tc[-1]
        return c

    def build(self, keys):
        # Initialize
        self.tws = tws
        self.head = None
        self.tail = None
        # Build by keys
        for key in keys:
            try:
                self.append(key)
            except:
                self.head = None
                self.tail = None
                raise KeyNotExistError("%s does not exist" % key)
        return

    def traverse(self):
        traverse = []
        n = self.head
        if (self.head == None):
            return []
        while (n.succ != None):
            traverse.append(n)
            n = n.succ
        traverse.append(n)
        return traverse

    def append(self, key):
        # NOTE: Does not need to update time windows
        n = JobNode(key)
        # interval will not be None since the last interval is to inf
        # self.tws[key]['twTree'].print()
        (t, interval) = self.tws[key]['twTree'].earlisetAvail(self.makespan)  # O(log I)
        # print(t, interval.value, key, self.makespan)
        n.ts = t
        n.te = t + interval.value
        n.prev = self.tail
        n.succ = None

        # For empty linked list
        if (self.head == None):
            self.head = n
            self.tail = n
        # If there is only one element
        elif (self.head != None and self.head.succ == None):
            self.head.succ = n
            self.tail = n
        else:
            self.tail.succ = n
            self.tail = n
        return

    def remove(self, key):
        if (self.head == None):
            raise KeyNotExistError("%s does not exist" % key)
        elif (self.head != None and self.head.succ == None):
            if (self.head.key == key):
                self.head = None
                self.tail = None
            else:
                raise KeyNotExistError("%s does not exist" % key)
        # At lease there are more than one node
        elif (self.head != None and self.head.succ != None):
            # Search from the beginning
            n = self.query(key)
            # Remove key node that is not the head
            if (n.prev != None):
                # Remove n
                n.prev.succ = n.succ
            else:
                self.head = n.succ
            if (n.succ != None):
                n.succ.prev = n.prev
                # NOTE: 移除了一个节点之后原先后续的节点开始要更新time windows, 直至第一个没有被更新的job
                self.updateFromNode(n.succ)
            else:
                self.tail = n.prev           
        return

    def insert(self, key, newKey):
        # insert newKey right BEFORE key, if key == None, append to the linked list
        # Old:  = = =  p   n s = =
        # New:  = = p newN n s = =
        if (self.head == None):
            raise KeyNotExistError("%s does not exist" % key)
        if (key == None):
            self.append(newKey)
        elif (self.head != None and self.head.succ == None):
            if (self.head.key == key):
                n = self.head
                newN = JobNode(newKey)
                (t, interval) = self.tws[newKey]['twTree'].earlisetAvail(0)
                newN.ts = t
                newN.te = t + interval.value
                newN.succ = n
                newN.prev = None
                self.head = newN
                self.tail = n
                n.prev = newN
                n.succ = None
                self.updateFromNode(n)
            else:
                raise KeyNotExistError("%s does not exist" % key)
        elif (self.head != None and self.head.succ != None):
            # Find key
            n = self.query(key)
            # Create a newKey node
            newN = JobNode(newKey)
            (t, interval) = self.tws[newKey]['twTree'].earlisetAvail(n.prev.te if n.prev != None else 0)  # O(log I)
            newN.ts = t
            newN.te = t + interval.value
            # Insert newN between n and n.succ
            newN.succ = n
            newN.prev = n.prev
            # Update the link
            if (n.prev != None):
                n.prev.succ = newN
            else:
                self.head = newN
            n.prev = newN
            self.updateFromNode(n)
        return

    def cheapInsert(self, key):
        if (self.head == None):
            self.append(key)
            return
        
        # Naive way
        # Step 1: Insert to head of linked list
        # Step 2: Swap from head to tail, log the makespan
        # Step 3: Remove the tail
        # Step 4: Insert to best position

        # Step 1: Insert to head of linked list
        bestInsert = self.head.key
        self.insert(self.head.key, key)
        bestMakespan = self.makespan
        # print("Step 1")
        # self.print()

        # Step 2: Swap from head to tail, log the makespan
        n = self.head.succ # Must be not-None, the origin head
        pos = n.key
        while(n != None):
            insert = n.succ.key if n.succ != None else None
            self.swap(key, pos)
            # print("Step 2", pos)
            # self.print()
            if (self.makespan < bestMakespan):
                bestInsert = insert
                bestMakespan = self.makespan
            n = n.succ.succ
            if (n != None):
                pos = n.key

        if (bestInsert != None):
            # Step 3: Remove the tail
            self.tail.prev.succ = None
            self.tail = self.tail.prev
            # print("Step 3")
            # self.print()

            # Step 4: Insert to the best position
            self.insert(bestInsert, key)
            # print("Step 4")
            # self.print()
        else:
            # If `bestInsert == None` we don need to remove the tail and add it back
            pass
        return

    def replace(self, oldKey, newKey):
        if (self.head == None):
            raise KeyNotExistError("%s does not exist" % key)
        elif (self.head != None and self.head.succ == None):
            if (self.head.key == oldKey):
                n = JobNode(newKey)
                # interval will not be None since the last interval is to inf
                (t, interval) = self.tws[newKey]['twTree'].earlisetAvail(0)  # O(log I)
                n.ts = t
                n.te = t + interval.value
                n.prev = None
                n.succ = None
                self.head = n
                self.tail = n
            else:
                raise KeyNotExistError("%s does not exist" % oldKey)
        # At lease there are more than one node
        elif (self.head != None and self.head.succ != None):
            # Find key
            n = self.query(oldKey)
            # Create a newKey node
            newN = JobNode(newKey)
            (t, interval) = self.tws[newKey]['twTree'].earlisetAvail(n.prev.te if n.prev != None else 0)  # O(log I)
            newN.ts = t
            newN.te = t + interval.value
            # Replace
            newN.prev = n.prev
            newN.succ = n.succ
            # Update link
            if (n.prev != None):
                n.prev.succ = newN
            else:
                self.head = newN
            if (n.succ != None):
                n.succ.prev = newN
            else:
                self.tail = newN
            if (newN.succ != None):
                # NOTE: newN的时间窗已经计算过了，所以更新时间窗从newN的下一个开始
                self.updateFromNode(newN.succ)
        return

    def swap(self, keyI, keyJ):
        if (self.head == None):
            raise KeyNotExistError("%s and %s does not exist." % (keyI, keyJ))
        # Find nI and nJ
        # FIXME: these two lines can be merged into one query
        nI = None
        nJ = None
        n = self.head
        while (n != None):
            if (n.key == keyI):
                nI = n
            if (n.key == keyJ):
                nJ = n
            if (nI != None and nJ != None):
                break
            n = n.succ

        if (n == None):
            if (nI == None and nJ != None):
                raise KeyNotExistError("%s does not exist" % keyI)
            if (nI != None and nJ == None):
                raise KeyNotExistError("%s does not exist" % keyJ)
            if (nI == None and nJ == None):
                raise KeyNotExistError("%s and %s does not exist" % (keyI, keyJ))

        if ((nI.succ != None and nI.succ.key == keyJ)
            or (nJ.succ != None and nJ.succ.key == keyI)):
            # Old: = = = p nI nJ s = = =
            # New: = = = p nJ nI s = = =
            if (nI.succ == nJ):
                if (nI.prev != None):
                    nI.prev.succ = nJ
                else:
                    self.head = nJ
                if (nJ.succ != None):
                    nJ.succ.prev = nI
                else:
                    self.tail = nI

                nI.prev, nI.succ, nJ.prev, nJ.succ = nJ, nJ.succ, nI.prev, nI
            else:
                if (nJ.prev != None):
                    nJ.prev.succ = nI
                else:
                    self.head = nI
                if (nI.succ != None):
                    nI.succ.prev = nJ
                else:
                    self.tail = nJ
                nJ.prev, nJ.succ, nI.prev, nI.succ = nI, nI.succ, nJ.prev, nJ
        else:
            # Update links
            # Old: = = = pI nI sI = = pJ nJ sJ = = =
            # New: = = = pI nJ sI = = pJ nI sJ = = =
            if (nI.prev != None):
                nI.prev.succ = nJ
            else:
                self.head = nJ
            if (nI.succ != None):
                nI.succ.prev = nJ
            else:
                self.tail = nJ
            if (nJ.prev != None):
                nJ.prev.succ = nI
            else:
                self.head = nI
            if (nJ.succ != None):
                nJ.succ.prev = nI
            nI.prev, nJ.prev = nJ.prev, nI.prev
            nI.succ, nJ.succ = nJ.succ, nI.succ
        
        # Update links
        # FIXME: this needs to be improved
        self.updateFromNode(nJ)
        if (nJ.succ != None):
            self.updateFromNode(nJ.succ)
        self.updateFromNode(nI)
        if (nI.succ != None):
            self.updateFromNode(nI.succ)
        return

    def updateFromNode(self, n=None):
        # Update the time windows starts from n, until we find a tw that is invariant
        if (n == None):
            n = self.head
        while (n != None):
            (t, interval) = self.tws[n.key]['twTree'].earlisetAvail(n.prev.te if n.prev != None else 0)
            newTs = t
            newTe = t + interval.value
            if (abs(newTs - n.ts) < CONST_EPSILON and abs(newTe - n.te) < CONST_EPSILON):
                break
            else:
                n.ts = newTs
                n.te = newTe
                n = n.succ
        return

