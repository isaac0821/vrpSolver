import math
from .common import *

class RingNode(object):
    def __init__(self, key, value=None, prev: 'RingNode' = None, next: 'RingNode' = None):
        self.key = key
        self.value = value if value != None else key
        self.prev = prev if prev != None else RingNilNode()
        self.next = next if next != None else RingNilNode()

    def clone(self):
        newObj = RingNode(self.key, self.value)
        newObj.prev = self.prev
        newObj.next = self.next
        return newObj

    @property
    def isNil(self):
        return False

    def __repr__(self):
        s =("{key: " + str(self.key) + ", "
            + "value: " + str(self.value) + ", "
            + "prev: " + (str(self.prev.key) if (not self.prev.isNil) else "None") + ", "
            + "next: " + (str(self.next.key) if (not self.next.isNil) else "None") + "} ")
        return s

class RingNilNode(RingNode):
    def __init__(self):
        return

    @property
    def isNil(self):
        return True

class Ring(object):
    def __init__(self):
        self.head = RingNilNode()
        self._count = 0

    @property
    def isEmpty(self):
        if (self.head.isNil):
            return True
        else:
            return False

    @property
    def count(self):
        return self._count

    def __repr__(self):
        return self.traverse().__repr__()

    def rehead(self, key):
        self.head = self.query(key)

    def query(self, key) -> "RingNode":
        if (self.head.isNil):
            raise EmptyError("ERROR: The route is empty.")
        cur = self.head
        trvFlag = False
        queried = 0
        while (not trvFlag):
            if (cur.key == key):
                return cur
            else:
                cur = cur.next
                queried += 1
                if (queried > self._count):
                    raise OutOfRangeError("ERROR: Unexpected loop")
                if (cur.key == self.head.key):
                    trvFlag = True
        if (trvFlag):
            return RingNilNode()

    def traverse(self, closeFlag=False) -> list:
        route = []
        cur = self.head
        while (not cur.isNil):
            if (len(route) > self._count):
                raise OutOfRangeError("ERROR: Unexpected loop")
            route.append(cur)
            cur = cur.next
            if (cur.key == self.head.key):
                break
        if (closeFlag):
            route.append(self.head)
        return route

    def insert(self, m, n):
        if (n.isNil):
            raise EmptyError("ERROR: Cannot insert an empty node.")
        if (self.head.isNil):
            self.head = n
            n.next = n
            n.prev = n
            self._count += 1
        else:
            m.next.prev = n
            n.prev = m
            n.next = m.next
            m.next = n
            self._count += 1

    def append(self, n):
        if (n.isNil):
            raise EmptyError("ERROR: Cannot insert an empty node.")
        if (self.head.isNil):
            self.head = n
            n.next = n
            n.prev = n
            self._count += 1
        else:
            return self.insert(self.head.prev, n)

    def remove(self, n):
        n.prev.next = n.next
        n.next.prev = n.prev
        self._count -= 1

class RouteNode(RingNode):
    def __init__(self, key, value=None, prev: 'RouteNode'=None, next: 'RouteNode'=None):
        self.key = key
        self.value = value if value != None else key
        self.prev = prev if prev != None else RouteNilNode()
        self.next = next if next != None else RouteNilNode()

class RouteNilNode(RouteNode):
    def __init__(self):
        return

    @property
    def isNil(self):
        return True

class Route(Ring):
    def __init__(self, tau, asymFlag=False):
        self.head = RouteNilNode()
        self.tau = tau
        self.dist = 0
        self._revDist = 0
        self.asymFlag = asymFlag
        self._count = 0

    @property
    def revDist(self):
        if (self.asymFlag):
            return self._revDist
        else:
            return self.dist

    def reverse(self):
        tra = self.traverse()
        for n in tra:
            n.prev, n.next = n.next, n.prev
        if (self.asymFlag):
            self.dist, self._revDist = self._revDist, self.dist

    def rotate(self, s, e):
        # = = s.prev s s.next = = = e.prev e e.next = = 
        # = = s.prev e e.prev = = = s.next s e.next = = 

        sPrev = s.prev
        eNext = e.next

        s.prev.next = e
        e.next.prev = s

        # 如果是asymmetric
        distS2E = 0
        distE2S = 0
        # 计算 s => e之间累计的距离
        if (self.asymFlag):
            k = s
            while (k.key != e.key):
                distS2E += self.tau[k.key, k.next.key]
                distE2S += self.tau[k.next.key, k.key]
                k = k.next

        tra = []
        k = s
        while (k.key != e.key):
            tra.append(k)
            k = k.next
        tra.append(e)

        for n in tra:
            n.prev, n.next = n.next, n.prev

        e.prev = sPrev
        s.next = eNext

        self.dist = self.dist - self.tau[sPrev.key, s.key] - self.tau[e.key, eNext.key] + self.tau[sPrev.key, e.key] + self.tau[s.key, eNext.key] - distS2E + distE2S
        self._revDist = self._revDist - self.tau[eNext.key, e.key] - self.tau[s.key, sPrev.key] + self.tau[eNext.key, s.key] + self.tau[e.key, sPrev.key] - distE2S + distS2E

    def query(self, key) -> "RouteNode":
        if (self.head.isNil):
            raise EmptyError("ERROR: The route is empty.")
        cur = self.head
        trvFlag = False
        while (not trvFlag):
            if (cur.key == key):
                return cur
            else:
                cur = cur.next
                if (cur.key == self.head.key):
                    trvFlag = True
        if (trvFlag):
            return RouteNilNode()

    def insert(self, m, n):
        if (n.isNil):
            raise EmptyError("ERROR: Cannot insert an empty node.")
        if (self.head.isNil):
            self.head = n
            n.next = n
            n.prev = n
            self.dist = 0
            if (self.asymFlag):
                self._revDist = 0
            self._count = 1
        else:
            self.dist += self.tau[m.key, n.key] + self.tau[n.key, m.next.key] - self.tau[m.key, m.next.key]
            if (self.asymFlag):
                self._revDist += self.tau[m.next.key, n.key] + self.tau[n.key, m.key] - self.tau[m.next.key, m.key]
            m.next.prev = n
            n.prev = m
            n.next = m.next
            m.next = n            
            self._count += 1
        return

    def append(self, n):
        if (self.head.isNil):
            self.head = n
            n.next = n
            n.prev = n
            self.dist = 0
            if (self.asymFlag):
                self._revDist = 0
            self._count = 1
        else:
            return self.insert(self.head.prev, n)

    def remove(self, n):
        self.dist += self.tau[n.prev.key, n.next.key] - self.tau[n.prev.key, n.key] - self.tau[n.key, n.next.key]
        if (self.asymFlag):
            self._revDist += self.tau[n.next.key, n.prev.key] - self.tau[n.key, n.prev.key] - self.tau[n.next.key, n.key]
        n.prev.next = n.next
        n.next.prev = n.prev
        if (self.head.key == n.key):
            self.head = n.next
        self._count -= 1

    def swap(self, n):
        nPrev = n.prev
        nNext = n.next
        nNNext = n.next.next
        # Pointers
        nPrev.next = nNext
        n.prev = nNext
        n.next = nNNext
        nNext.prev = nPrev
        nNext.next = n
        nNNext.prev = n
        # Calculate dist
        self.dist += (self.tau[nPrev.key, nNext.key] + self.tau[nNext.key, n.key] + self.tau[n.key, nNNext.key]
                    - self.tau[nPrev.key, n.key] - self.tau[n.key, nNext.key] - self.tau[nNext.key, nNNext.key])
        if (self.asymFlag):
            self._revDist += (self.tau[nNNext.key, n.key] + self.tau[n.key, nNext.key] + self.tau[nNext.key, nPrev.key]
                            - self.tau[nNNext.key, nNext.key] + self.tau[nNext.key, n.key] - self.tau[n.key, nPrev.key])

    def exchange(self, nI, nJ):
        if (nI == nJ):
            raise KeyExistError("ERROR: Cannot swap itself")
        # Old: = = = i j k l m n = = = | -> exchange(j, m)
        # New: = = = i m k l j n = = =
        if (nI.next.key == nJ.key):
            self.swap(nI)
        if (nI.next.next.key == nJ.key):
            # Old: = = pI nI nX nJ sJ = =
            # New: = = pI nJ nX nI sJ = =
            pI = nI.prev
            nX = nI.next
            sJ = nJ.next
            pI.next = nJ
            nJ.prev = pI
            nJ.next = nX
            nX.prev = nJ
            nX.next = nI
            nI.prev = nX
            nI.next = sJ
            sJ.prev = nI
            self.dist += (self.tau[pI.key, nJ.key] + self.tau[nJ.key, nX.key] + self.tau[nX.key, nI.key] + self.tau[nI.key, sJ.key]
                - self.tau[pI.key, nI.key] - self.tau[nI.key, nX.key] - self.tau[nX.key, nJ.key] - self.tau[nJ.key, sJ.key])
            if (self.asymFlag):
                self._revDist += (self.tau[sJ.key, nI.key] + self.tau[nI.key, nX.key] + self.tau[nX.key, nJ.key] + self.tau[nJ.key, pI.key]
                    - self.tau[sJ.key, nJ.key] - self.tau[nJ.key, nX.key] - self.tau[nX.key, nI.key] - self.tau[nI.key, pI.key])

        else:
            # Old: = = pI nI sI x x x pJ nJ sJ = =
            # New: = = pI nJ sI x x x pJ nI sJ = =
            pI = nI.prev
            sI = nI.next
            pJ = nJ.prev
            sJ = nJ.next            
            pI.next = nJ
            nJ.prev = pI
            nJ.next = sI
            sI.prev = nJ
            pJ.next = nI
            nI.prev = pJ
            nI.next = sJ
            sJ.prev = nI
            self.dist += (self.tau[pI.key, nJ.key] + self.tau[nJ.key, sI.key] + self.tau[pJ.key, nI.key] + self.tau[nI.key, sJ.key]
                - self.tau[pI.key, nI.key] - self.tau[nI.key, sI.key] - self.tau[pJ.key, nJ.key] - self.tau[nJ.key, sJ.key])
            if (self.asymFlag):
                self._revDist += (self.tau[sJ.key, nI.key] + self.tau[nI.key, pJ.key] + self.tau[sI.key, nJ.key] + self.tau[nJ.key, pI.key]
                    - self.tau[sJ.key, nJ.key] - self.tau[nJ.key, pJ.key] - self.tau[sI.key, nI.key] - self.tau[nI.key, pI.key])

    def cheapestInsert(self, n):
        # First, insert it after self.head
        # Before: ... --> head -> xxx -> head.next --> ...
        # After:  ... --> head ->  n  -> head.next --> ...
        self.insert(self.head, n)
        sofarCheapestCost = self.dist if not self.asymFlag else min(self.dist, self._revDist)
        sofarCheapestKey = self.head
        if (self._count <= 2):
            return
        cur = self.head
        trvFlag = True
        while (trvFlag):
            self.swap(n)
            cur = cur.next
            newCost = self.dist if not self.asymFlag else min(self.dist, self._revDist)
            if (newCost < sofarCheapestCost):
                sofarCheapestCost = newCost
                sofarCheapestKey = cur
            if (cur.key == self.head.key):
                trvFlag = False
        self.remove(n)
        self.insert(sofarCheapestKey, n)
        if (self.asymFlag and self.dist > self._revDist):
            self.reverse()

    def findLargestRemoval(self, noRemoval=None) -> int:
        if (noRemoval == None):
            noRemoval = [self.head.key]
        bestRemovalCost = self.dist if not self.asymFlag else max(self.dist, self._revDist)
        bestRemovalKey = None
        cur = self.head
        trvFlag = True
        while (trvFlag):
            cur = cur.next
            if (cur.key not in noRemoval):
                newDist = (self.dist + self.tau[cur.prev.key, cur.next.key] 
                    - self.tau[cur.prev.key, cur.key] - self.tau[cur.key, cur.next.key])
                if (self.asymFlag):
                    newRevDist = (self._revDist + self.tau[cur.next.key, cur.prev.key]
                        - self.tau[cur.next.key, cur.key] - self.tau[cur.key, cur.prev.key])
                    newDist = min(newDist, newRevDist)
                if (newDist < bestRemovalCost):
                    bestRemovalCost = newDist
                    bestRemovalKey = cur.key
            if (cur.key == self.head.key):
                trvFlag = False
        return {
            'bestRemovalKey': bestRemovalKey,
            'bestRemovalCost': bestRemovalCost
        }

    def impv2Opt(self):
        # NOTE: Now its better.
        oriHeadKey = self.head.key
        nI = self.head.next
        sofarBestDist = self.dist if not self.asymFlag else min(self.dist, self._revDist)
        improvedFlag = False
        canImpvFlag = True
        while (canImpvFlag):
            canImpvFlag = False

            endKey = nI.prev.key
            while (nI.key != endKey):
                # 1. First attempt
                # Old: = = = nI nINext nJ nJNext nJ2Next nJ3Next nJ4Next = = = | -> exchange(nINext, nJ)        | -> exchange(nINext, nJ)
                # 2. Follow-up attempts
                # New: = = = nI nJ nINext nJNext nJ2Next nJ3Next nJ4Next = = = | -> exchange(nINext, nJNext)    | -> remove(nJNext)
                #                                                                -> exchange(nJ, nJNext)          -> insert(nI, nJNext)
                # New: = = = nI nJNext nJ nINext nJ2Next nJ3Next nJ4Next = = = | -> exchange(nINext, n2JNext)   | -> remove(nJ2Next)
                #                                                                -> exchange(nJ, nJ2Next)         -> insert(nI, nJ2Next)
                #                                                                -> exchange(nJNext, nJ2Next)
                # New: = = = nI nJ2Next nJNext nJ nINext nJ3Next nJ4Next = = = | -> exchange(nINext, n3JNext)   | -> remove(nJ3Next)
                #                                                                -> exchange(nJ, nJ3Next)         -> insert(nI, nJ3Next)
                #                                                                -> exchange(nJNext, nJ3Next) 
                #                                                                -> exchange(nJ2Next, nJ3Next)
                # New: = = = nI nJ3Next nJ2Next nJNext nJ nINext nJ4Next = = =
                # ...
                # 3. Recover to initial status (Until: nJXNext == nIPrev)
                # Old: nI nJ = = = nINext nJNext (nIPrev) | -> exchange(nJNext, nI)
                # New: nI nJNext nJ = = = nINext
                # self.reverse()

                # 1. First attempt
                nJ = nI.next.next
                nINext = nI.next
                self.swap(nINext)
                newDist = self.dist if not self.asymFlag else min(self.dist, self._revDist)
                if (newDist < sofarBestDist):
                    sofarBestDist = newDist
                    self.rehead(oriHeadKey)
                    canImpvFlag = True
                    improvedFlag = True
                    break

                # 2. Follow-up attempts
                for _ in range(self._count - 4):
                    nJXNext = nINext.next
                    self.remove(nJXNext)
                    self.insert(nI, nJXNext)
                    newDist = self.dist if not self.asymFlag else min(self.dist, self._revDist)
                    if (newDist < sofarBestDist):
                        sofarBestDist = newDist
                        self.rehead(oriHeadKey)
                        canImpvFlag = True
                        improvedFlag = True
                        break
                if (canImpvFlag):
                    break

                # 3. Recover to initial status
                nJXNext = nINext.next
                self.swap(nJXNext)
                self.reverse()
                self.rehead(oriHeadKey)
                
                nI = nI.next
        return improvedFlag

class TreeNode(object):
    def __init__(self, key, value, parent: 'TreeNode' = None, treeNodes: list['TreeNode'] = None):
        self.key = key
        self.value = value
        self.parent = parent if parent != None else TreeNilNode()
        self.treeNodes = treeNodes if treeNodes != None else [TreeNilNode()]

    @property
    def isNil(self):
        return False

    @property
    def isChildren(self):
        if (len(self.treeNodes) == 1 and self.treeNodes[0].isNil):
            return True
        else:
            return False

class TreeNilNode(TreeNode):
    def __init__(self):
        return

    @property
    def isNil(self):
        return True

class Tree(object):
    def __init__(self):
        self.nil = TreeNilNode()
        self.root = self.nil

    def __repr__(self):
        tr = self.traverse()
        return str(tr)

    @property
    def isEmpty(self):
        return self.root.isNil

    def traverse(self):
        if (self.root.isNil):
            return []
        else:
            tra = [self.root]
        tra.extend(self._traverseBreath(self.root))
        return tra

    def _traverseBreath(self, n):
        tra = []
        for treeNode in n.treeNodes:
            if (not treeNode.isNil):
                tra.append(treeNode)
        for treeNode in n.treeNodes:
            if (not treeNode.isNil):
                tra.extend(self._traverseBreath(treeNode))
        return tra

    # Query using key
    def query(self, key):
        searched = self._search(self.root, key)
        if (searched != None):
            return searched
        else:
            raise KeyNotExistError("ERROR: %s does not exist." % key)
    
    def _search(self, n, key):
        if (n.isNil):
            return None
        else:
            if (n.key == key):
                return n
            for treeNode in n.treeNodes:
                searched = self._search(treeNode, key)
                if (searched != None and not searched.isNil):
                    return searched
        
    def insert(self, n, treeNode = None):
        if (treeNode == None):
            if (self.root.isNil):
                self.root = n
                return
            else:
                raise KeyExistError("ERROR: Root exists.")

        if (n.isChildren):
            n.treeNodes = [treeNode]
        else:
            n.treeNodes.append(treeNode)
        treeNode.parent = n
        return

    def traverseChildren(self):
        if (self.root.isNil):
            return []
        else:
            if (self.root.isChildren):
                return [self.root]
        children = self._traverseChildren(self.root)
        return children

    def _traverseChildren(self, n):
        tra = []
        for treeNode in n.treeNodes:
            if (treeNode.isChildren):
                tra.append(treeNode)
        for treeNode in n.treeNodes:
            if (not treeNode.isChildren):
                tra.extend(self._traverseChildren(treeNode))
        return tra

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
    def isEmpty(self):
        return self.root.isNil

    @property
    def count(self):
        return len(self.traverse())

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

class SimpleTimeWindowsTree(RedBlackTree):
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
    
class IntervalTree(RedBlackTree):
    def __init__(self):
        self.nil = IntervalTreeNilNode()
        self.root = self.nil

    def querySingular(self, t) -> list[IntervalTreeNode]:
        return

    def queryInterval(self, ts, te) -> list[IntervalTreeNode]:
        return

    def earliestNextInterval(self, t) -> IntervalTreeNode:
        return

    def earliestNextAvailable(self, t) -> float:
        return

    def _search(self, n, t):
        if (n == self.nil):
            return self.nil
        elif (n.lower <= t and t <= n.upper):
            return n
        elif (t < n.lower):
            return self._search(n.left, t)
        elif (t > n.upper):
            return self._search(n.right, t)

class Job(object):
    def __init__(self, length):
        self.length = length

class JobTimeWindowed(Job):
    def __init__(self, key, twTree):
        self.key = key
        self.twTree = twTree
        
class JobNode(object):
    def __init__(self, key, job, value=None, ts=None, te=None, prev=None, succ=None):
        self.key = key
        self.job = job
        self.value = value
        self.ts = ts
        self.te = te
        self.prev = prev
        self.succ = succ

class ScheduleList(object):
    def __init__(self, tabuTWs = None):
        self.head = None
        self.tail = None         # tail要单独处理, 只有一个元素的时候tail是head的复制
        self.tabuTWs = tabuTWs   # 禁止时间窗

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
        c = ScheduleList(self.tabuTWs)
        tr = self.traverse()
        tc = []
        for i in range(len(tr)):
            tc.append(JobNode(
                key = tr[i].key,
                value = tr[i].value,
                ts = tr[i].ts,
                te = tr[i].te,
                prev = tc[i - 1] if i >= 1 else None))
        for i in range(len(tr)):
            tc[i].succ = tc[i + 1] if (i < len(tr) - 1) else None
        if (len(tc) > 0):
            c.head = tc[0]
            c.tail = tc[-1]
        return c

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

    def traverseKey(self):
        traverse = []
        n = self.head
        if (self.head == None):
            return []
        while (n.succ != None):
            traverse.append(n.key)
            n = n.succ
        traverse.append(n.key)
        return traverse

    def append(self, key, job):
        # NOTE: Does not need to update time windows
        n = JobNode(key, job)
        n.ts = self.te
        n.te = self.te + n.job.length
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

    def insert(self, key, newJID):
        # insert newJID right BEFORE key, if key == None, append to the linked list
        # Old:  = = =  p   n s = =
        # New:  = = p newN n s = =
        if (self.head == None):
            raise KeyNotExistError("%s does not exist" % key)
        if (key == None):
            self.append(newJID)
        elif (self.head != None and self.head.succ == None):
            if (self.head.key == key):
                n = self.head
                newN = JobNode(newJID)
                (t, interval) = self.jobs[newJID]['twTree'].earlisetAvail(0)
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
            # Create a newJID node
            newN = JobNode(newJID)
            (t, interval) = self.jobs[newJID]['twTree'].earlisetAvail(n.prev.te if n.prev != None else 0)  # O(log I)
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
            self.exchange(key, pos)
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

    def replace(self, oldJID, newJID):
        if (self.head == None):
            raise KeyNotExistError("%s does not exist" % key)
        elif (self.head != None and self.head.succ == None):
            if (self.head.key == oldJID):
                n = JobNode(newJID)
                # interval will not be None since the last interval is to inf
                (t, interval) = self.jobs[newJID]['twTree'].earlisetAvail(0)  # O(log I)
                n.ts = t
                n.te = t + interval.value
                n.prev = None
                n.succ = None
                self.head = n
                self.tail = n
            else:
                raise KeyNotExistError("%s does not exist" % oldJID)
        # At lease there are more than one node
        elif (self.head != None and self.head.succ != None):
            # Find key
            n = self.query(oldJID)
            # Create a newJID node
            newN = JobNode(newJID)
            (t, interval) = self.jobs[newJID]['twTree'].earlisetAvail(n.prev.te if n.prev != None else 0)  # O(log I)
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

    def exchange(self, keyI, keyJ):
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
            (t, interval) = self.jobs[n.key]['twTree'].earlisetAvail(n.prev.te if n.prev != None else 0)
            newTs = t
            newTe = t + interval.value
            if (abs(newTs - n.ts) < CONST_EPSILON and abs(newTe - n.te) < CONST_EPSILON):
                break
            else:
                n.ts = newTs
                n.te = newTe
                n = n.succ
        return

