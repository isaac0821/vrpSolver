from .common import *
from .error import *

class BSTreeNode:
    def __init__(self, key, value, parent=None, left=None, right=None):
        self.key = key
        self.value = value
        self.parent = parent
        self.left = left
        self.right = right

class BSTree:
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
    def __init__(self, lower, upper, value, parent=None, left=None, right=None):
        self.key = lower
        self.value = value
        self.lower = lower # float('-inf') represents infinite
        self.upper = upper # float('inf') represents infinite
        self.parent = parent
        self.left = left
        self.right = right
        self.color = None  # RedBlackTree

class IntervalTree(RedBlackTree):
    # Interval tree ===========================================================
    # NOTE: Interval tree is an augmented version of Red Black tree
    # NOTE: Used in time windows calculation
    def __init__(self):
        self.nil = RedBlackTreeNode(key=None, value=None, color='B')
        self.root = self.nil

    def query(self, t):
        intervals = []
        def search(n, t):
            if (n == self.nil):
                return self.nil
            if (n.lower <= t and t <= n.upper):
                intervals.append(n)
            if (n.key == t):
                return n
            if (t < n.key):
                return search(n.left, t)
            else:
                return search(n.right, t)
        search(self.root, t)
        return intervals

    def earlisetAvail(self, t):


        n = self.query(key)
        if (n != self.nil):
            return self.next(n)


        if (n.right != self.nil):
            return self.min(n.right)
        else:
            y = n.parent
            while (y != self.nil and n == y.right):
                n = y
                y = y.parent
            return y

        return nextAvail

class LinkedListNode:
    def __init__(self, key, value, prev=None, succ=None):
        self.key = key
        self.value = value
        self.prev = prev
        self.succ = succ

class LinkedList:
    # Linked list =============================================================
    # NOTE: used in TSP/VRP, representing routes
    def __init__(self):
        self.nil = None
        self.head = self.nil
        self.tail = self.nil

    def prepend(self, n):
        n.succ = self.head
        n.prev = self.nil
        if (self.head != self.nil):
            self.head.prev = n
        self.head = n
        return

    def append(self, n):
        if (self.head == self.nil):
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
        if (x.succ != self.nil):
            x.succ.prev = n
        x.succ = n
        return

    def deleteByKey(self, key):
        n = self.query(key)
        self.delete(n)
        return
    def delete(self, n):
        if (n.prev != self.nil):
            n.prev.succ = n.succ
        else:
            self.head = n.succ
        if (n.succ != self.nil):
            n.succ.prev = n.prev
        return

    def query(self, key):
        n = self.head
        while (n != self.nil and n.key != key):
            n = n.succ
        return n

    def traverse(self, mode="Left"):
        traverse = []
        if (mode == "Left"):
            n = self.head
            while (n.succ != self.nil):
                traverse.append(n)
                n = n.succ
            traverse.append(n)
        elif (mode == "Right"):
            n = self.tail
            while (n.prev != self.nil):
                traverse.append(n)
                n = n.prev
            traverse.append(n)
        return traverse
