import math
from .error import *

class Ring():
    def _init__(self, l:list=[]):
        self.element = [i for i in l]

    def append(self, ele):
        self.element.append(ele)

    def between(self, s, e):
        btw = []
        if (s not in self.element):
            raise KeyNotExistError(s)
        if (e not in self.element):
            raise KeyNotExistError(e)
        idxS = self.element.index(s)
        idxE = self.element.index(e)

        if (idxS < idxE):
            for i in range(idxS, idxE):
                btw.append(self.element[i])
        else:
            for i in range(idxS, len(self.element)):
                btw.append(self.element[i])
            for i in range(idxE):
                btw.append(self.element[i])
        return btw

    def remove(self, ele):
        self.element.remove(ele)

    def next(self, ele):
        if (ele not in self.element):
            raise KeyNotExistError(ele)
        idx = self.element.index(ele)
        if (idx == len(self.element) - 1):
            return self.element[0]
        else:
            return self.element[idx + 1]

    def prev(self, ele):
        if (ele not in self.element):
            raise KeyNotExistError(ele)
        idx = self.element.index(ele)
        if (idx == 0):
            return self.element[-1]
        else:
            return self.element[idx - 1]

class BSTreeNode(object):
    def __init__(self, key:int, value, parent:'BSTreeNode'=None, left:'BSTreeNode'=None, right:'BSTreeNode'=None):
        self.key = key
        self.value = value
        self.parent = parent if parent != None else BSTreeNilNode()
        self.left = left if left != None else BSTreeNilNode()
        self.right = right if right != None else BSTreeNilNode()

    @property
    def isNil(self):
        return False

    def __repr__(self):
        s =("{key: " + str(self.key) + ", "
            + "value: " + str(self.value) + ", "
            + "parent: " + (str(self.parent.key) if (not self.parent.isNil) else "None") + ", "
            + "left: " + (str(self.left.key) if (not self.left.isNil) else "None") + ", "
            + "right: " + (str(self.right.key) if (not self.right.isNil) else "None") + "}\n")
        return s

class BSTreeNilNode(BSTreeNode):
    def __init__(self):
        return

    @property
    def isNil(self):
        return True

class BSTree(object):
    # Initialize with an empty tree
    def __init__(self):
        self.nil = BSTreeNilNode()
        self.root = self.nil


    def __repr__(self):
        tr = self.traverse()
        return str(tr)

    @property
    def count(self):
        return len(self.traverse())

    @property
    def isEmpty(self):
        return self.root == self.nil

    # Query using key
    def query(self, key:int):
        return self._search(self.root, key)
    def _search(self, n, key:int):
        if (n.isNil):
            return self.nil
        else:
            if (n.key == key):
                return n
            if (key < n.key):
                return self._search(n.left, key)
            else:
                return self._search(n.right, key)

    def traverse(self, mode='Left'):
        traverse = []
        def leftTraverse(n):
            if (not n.isNil):
                leftTraverse(n.left)
                traverse.append(n)
                leftTraverse(n.right)
            return
        def midTraverse(n):
            if (not n.isNil):
                traverse.append(n)
                midTraverse(n.left)                
                midTraverse(n.right)
            return
        def rightTraverse(n):
            if (not n.isNil):
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
    def traverseValue(self, mode='Left'):
        traverse = []
        def leftTraverse(n):
            if (not n.isNil):
                leftTraverse(n.left)
                traverse.append(n.value)
                leftTraverse(n.right)
            return
        def midTraverse(n):
            if (not n.isNil):
                traverse.append(n.value)
                midTraverse(n.left)                
                midTraverse(n.right)
            return
        def rightTraverse(n):
            if (not n.isNil):
                rightTraverse(n.right)
                traverse.append(n.value)
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
        if (n.isNil):
            raise EmptyError("ERROR: Cannot insert and Nil node to the tree")

        x = self.root
        y = self.nil
        while (not x.isNil):
            y = x
            if (n.key < x.key):
                x = x.left
            elif(n.key > x.key):
                x = x.right
            else:
                raise KeyExistError("ERROR: Key %s exists in BST." % n.key)

        n.parent = y
        if (y.isNil):
            self.root = n
        elif (n.key < y.key):
            y.left = n
        else:
            y.right = n
        return

    # Delete a BSTreeNode from BST
    def delete(self, key:int):
        n = self.query(key)
        if (not n.isNil):
            return self._delete(n)
        else:
            raise KeyNotExistError("ERROR: Cannot find key %s in BST" % key)
    def _delete(self, n:BSTreeNode):
        # Replace node u with node v in this previous location
        def _replace(u, v):
            if (u.parent.isNil):
                self.root = v
            elif (u == u.parent.left):
                u.parent.left = v
            else:
                u.parent.right = v
            if (not v.isNil):
                v.parent = u.parent

        if (n.left.isNil):
            _replace(n, n.right)
        elif (n.right.isNil):
            _replace(n, n.left)
        else:
            y = self.min(n.right)
            if (y != n.right):
                _replace(y, y.right)
                y.right = n.right
                y.right.parent = y
            _replace(n, y)
            y.left = n.left
            y.left.parent = y

        return

    # Predecessor
    def prev(self, key:int):
        n = self.query(key)
        if (n.isNil):
            return self._prev(n)
        else:
            raise KeyNotExistError("ERROR: Cannot find key %s in BST" % key)
    def _prev(self, n):
        if (not n.left.isNil):
            return self.max(n.left)
        else:
            y = n.parent
            while (not y.isNil and n == y.left):
                n = y
                y = y.parent
            return y

    # Successor
    def next(self, key:int):
        n = self.query(key)
        if (n.isNil):
            return self._next(n)
        else:
            raise KeyNotExistError("ERROR: Cannot find key %s in BST" % key)
    def _next(self, n):
        if (n.right.isNil):
            return self.min(n.right)
        else:
            y = n.parent
            while (not y.isNil and n == y.right):
                n = y
                y = y.parent
            return y

    # Get minimum
    def min(self, n):
        if (not n.isNil):
            while (not n.left.isNil):
                n = n.left
            return n
        else:
            raise EmptyError("ERROR: cannot perform min() on an empty tree.")

    # Get maximum
    def max(self, n):
        if (not n.isNil):
            while (not n.right.isNil):
                n = n.right
            return n
        else:
            raise EmptyError("ERROR: cannot perform min() on an empty tree.")

    # Left rotation
    def _leftRotation(self, x):
        y = x.right
        x.right = y.left
        if (not y.left.isNil):
            y.left.parent = x
        y.parent = x.parent
        if (x.parent.isNil):
            self.root = y
        elif (x == x.parent.left):
            x.parent.left = y
        else:
            x.parent.right = y
        y.left = x
        x.parent = y
        return

    # Right rotation
    def _rightRotation(self, x):
        y = x.left
        x.left = y.right
        if (not y.right.isNil):
            y.right.parent = x
        y.parent = x.parent
        if (x.parent.isNil):
            self.root = y
        elif (x == x.parent.right):
            x.parent.right = y
        else:
            x.parent.left = y
        y.right = x
        x.parent = y
        return

class RedBlackTreeNode(BSTreeNode):
    def __init__(self, key:int, value, parent:'RedBlackTreeNode'=None, left:'RedBlackTreeNode'=None, right:'RedBlackTreeNode'=None, color:str='B'):
        super().__init__(key, value)
        self.parent = parent if parent != None else RedBlackTreeNilNode()
        self.left = left if left != None else RedBlackTreeNilNode()
        self.right = right if right != None else RedBlackTreeNilNode()
        self.color = color

    def __repr__(self):
        s =("{key: " + str(self.key) + ", "
            + "value: " + str(self.value) + ", "
            + "color: " + str(self.color) + ", "
            + "parent: " + (str(self.parent.key) if (not self.parent.isNil) else "None") + ", "
            + "left: " + (str(self.left.key) if (not self.left.isNil) else "None") + ", "
            + "right: " + (str(self.right.key) if (not self.right.isNil) else "None") + "}")
        return s

class RedBlackTreeNilNode(RedBlackTreeNode):
    def __init__(self):
        self.color = 'B'

    @property
    def isNil(self):
        return True

class RedBlackTree(BSTree):
    def __init__(self):
        self.nil = RedBlackTreeNilNode()
        self.root = self.nil    

    # Insertion
    def insert(self, z):
        if (z.isNil):
            raise EmptyError("ERROR: cannot insert a Nil node to the tree.")

        x = self.root
        y = self.nil
        # Find parent of z
        while (not x.isNil):
            y = x
            if (z.key < x.key):
                x = x.left
            elif(z.key > x.key):
                x = x.right
            else:
                raise KeyExistError("Key %s exists in BST." % z.key)

        # Insert z to y
        z.parent = y
        if (y.isNil):
            self.root = z
        elif (z.key < y.key):
            y.left = z
        else:
            y.right = z
        z.left = self.nil
        z.right = self.nil
        z.color = 'R'
        while (z.parent.color == 'R'):
            if (z.parent == z.parent.parent.left):
                y = z.parent.parent.right
                if (y.color == 'R'):
                    z.parent.color ='B'
                    y.color = 'B'
                    z.parent.parent.color = 'R'
                    z = z.parent.parent
                else:
                    if (z == z.parent.right):
                        z = z.parent
                        super()._leftRotation(z)
                    z.parent.color = 'B'
                    z.parent.parent.color = 'R'
                    super()._rightRotation(z.parent.parent)
            else:
                y = z.parent.parent.left
                if (y.color == 'R'):
                    z.parent.color ='B'
                    y.color = 'B'
                    z.parent.parent.color = 'R'
                    z = z.parent.parent
                else:
                    if (z == z.parent.left):
                        z = z.parent
                        super()._rightRotation(z)
                    z.parent.color = 'B'
                    z.parent.parent.color = 'R'
                    super()._leftRotation(z.parent.parent)
        self.root.color = 'B'
        return

    def _delete(self, z):
        def _replace(u, v):
            if (u.parent.isNil):
                self.root = v
            elif (u == u.parent.left):
                u.parent.left = v
            else:
                u.parent.right = v
            v.parent = u.parent
        y = z
        yOriColor = y.color
        if (z.left.isNil):
            x = z.right
            _replace(z, z.right)
        elif (z.right.isNil):
            x = z.left
            _replace(z, z.left)
        else:
            y = super().min(z.right)
            yOriColor = y.color
            x = y.right
            if (y != z.right):
                _replace(y, y.right)
                y.right = z.right
                y.right.parent = y
            else:
                x.parent = y
            _replace(z, y)
            y.left = z.left
            y.left.parent = y
            y.color = z.color
        if (yOriColor == 'B'):
            while (x != self.root and x.color == 'B'):
                if (x == x.parent.left):
                    w = x.parent.right
                    if (w.color == 'R'):
                        w.color = 'B'
                        x.parent.color = 'R'
                        super()._leftRotation(x.parent)
                        w = x.parent.right
                    if (w.left.color == 'B' and w.right.color == 'B'):
                        w.color = 'R'
                        x = x.parent
                    else:
                        if (w.right.color == 'B'):
                            w.left.color = 'B'
                            w.color = 'R'
                            super()._rightRotation(w)
                            w = x.parent.right
                        w.color = x.parent.color
                        x.parent.color = 'B'
                        w.right.color = 'B'
                        super()._leftRotation(x.parent)
                        x = self.root
                else:
                    w = x.parent.left
                    if (w.color == 'R'):
                        w.color = 'B'
                        x.parent.color = 'R'
                        super()._rightRotation(x.parent)
                        w = x.parent.left
                    if (w.right.color == 'B' and w.left.color == 'B'):
                        w.color = 'R'
                        x = x.parent
                    else:
                        if (w.left.color == 'B'):
                            w.right.color = 'B'
                            w.color = 'R'
                            super()._leftRotation(w)
                            w = x.parent.left
                        w.color = x.parent.color
                        x.parent.color = 'B'
                        w.left.color = 'B'
                        super()._rightRotation(x.parent)
                        x = self.root
            x.color = 'B'
        return

class IntervalTreeNode(RedBlackTreeNode):
    def __init__(self, lower:float, upper:float, value, parent: 'IntervalTreeNode'=None, left: 'IntervalTreeNode' = None, right: 'IntervalTreeNode' = None):
        if (upper < lower):
            raise UnsupportedInputError("ERROR: `upper` should not be less than `lower`")
        super().__init__(key=1, value=value)
        self.key = (lower, upper)
        self.lower = lower
        self.upper = upper
        self.parent = parent if parent != None else IntervalTreeNilNode()
        self.left = left if left != None else IntervalTreeNilNode()
        self.right = right if right != None else IntervalTreeNilNode()

    def __repr__(self):
        s =("{lower: " + str(self.lower) + ", "
            + "upper: " + str(self.upper) + ", "
            + "value: " + str(self.value) + ", "
            + "parent: " + (str(self.parent.key) if (not self.parent.isNil) else "None") + ", "
            + "left: " + (str(self.left.key) if (not self.left.isNil) else "None") + ", "
            + "right: " + (str(self.right.key) if (not self.right.isNil) else "None") + "}")
        return s

    @property
    def max(self):
        return max(self.upper if self.upper != None else -float('inf'), 
            self.left.upper if self.left.upper != None else -float('inf'),
            self.right.upper if self.right.upper != None else -float('inf'))    

class IntervalTreeNilNode(IntervalTreeNode):
    def __init__(self):
        self.lower = None
        self.upper = None
        self.color = 'B'

    @property
    def isNil(self):
        return True
    
class IntervalTree(RedBlackTree):
    def __init__(self):
        self.nil = IntervalTreeNilNode()
        self.root = self.nil

    def querySingular(self, t) -> list[IntervalTreeNode]:

        return

    def queryInterval(self, ts, te) -> list[IntervalTreeNode]:

        return

    def earliestNext(self, t) -> float:
        en = float('inf')
        invs = self.querySingular(t)
        if (isinstance(invs, IntervalTreeNode)):
            en = t
        else:
            en = t

        return en

    @property
    def isSelfOverlapped(self):
        return
    def _checkOverlap(self, n) -> bool:
        if (n.left.isNil and n.right.isNil):
            return False

    def _search(self, n, t):
        if (n == self.nil):
            return self.nil
        elif (n.lower <= t and t <= n.upper):
            return n
        elif (t < n.lower):
            return self._search(n.left, t)
        elif (t > n.upper):
            return self._search(n.right, t)
