import math
from .error import *
from .const import *

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
        while (not trvFlag):
            if (cur.key == key):
                return cur
            else:
                cur = cur.next
                if (cur.key == self.head.key):
                    trvFlag = True
        if (trvFlag):
            return RingNilNode()

    def traverse(self, closeFlag=False) -> list:
        route = []
        cur = self.head
        while (not cur.isNil):
            route.append(cur)
            cur = cur.next
            if (cur.key == self.head.key):
                break
        if (closeFlag):
            route.append(self.head)
        return route

    def insert(self, key, n):
        if (n.isNil):
            raise EmptyError("ERROR: Cannot insert an empty node.")
        if (self.head.isNil):
            self.head = n
            n.next = n
            n.prev = n
            self._count += 1
        else:
            if (not self.query(n.key).isNil):
                raise KeyExistError("ERROR: %s already in the route." % n.key)
            m = self.query(key)
            if (m.isNil):
                raise KeyNotExistError("ERROR: %s does not exist." % key)
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
            return self.insert(self.head.prev.key, n)

    def remove(self, key):
        if (self.head.isNil):
            raise EmptyError("ERROR: The route is empty.")
        n = self.query(key)
        if (n.isNil):
            raise KeyNotExistError("ERROR: %s does not exist." % key)
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

    def reverseBetween(self, startKey, endKey):
        if (startKey == endKey):
            return
        tra = self.queryBetween(startKey, endKey)
        startPrev = tra[0].prev
        endNext = tra[-1].next
        # Calculate reverse distance
        s2eDist = 0
        e2sDist = 0
        if (self.asymFlag):
            cur = tra[0]
            while (cur != tra[-1]):
                s2eDist += self.tau[cur.key, cur.next.key]
                e2sDist += self.tau[cur.next.key, cur.key]
                cur = cur.next
        # Change pointers
        for n in tra:
            n.prev, n.next = n.next, n.prev
        tra[-1].prev = startPrev
        tra[0].next = endNext
        startPrev.next = tra[-1]
        endNext.prev = tra[0]
        # Update distances
        self.dist += (self.tau[startPrev.key, tra[-1].key] + self.tau[tra[0].key, endNext] 
            - self.tau[startPrev.key, tra[0].key] - self.tau[tau[-1].key, endNext.key])
        if (self.asymFlag):
            self._revDist += (self.tau[startPrev.key, tra[-1].key] + self.tau[tra[0].key, endNext] 
                - self.tau[startPrev.key, tra[0].key] - self.tau[tau[-1].key, endNext.key])
            self._revDist += e2sDist - s2eDist

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

    def queryBetween(self, startKey, endKey) -> list:
        if (self.head.isNil):
            raise EmptyError("ERROR: The route is empty.")
        cur = self.query(startKey)
        if (cur.isNil):
            return []
        trvFlag = False
        q = [cur]
        while (not trvFlag):
            if (cur.key == endKey):
                return q
            else:
                cur = cur.next
                q.append(cur)
                if (cur.key == endKey):
                    trvFlag = True
        return q

    def insert(self, key, n):
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
            if (not self.query(n.key).isNil):
                raise KeyExistError("ERROR: %s already in the route." % n.key)
            m = self.query(key)
            if (m.isNil):
                raise KeyNotExistError("ERROR: %s does not exist." % key)
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
            return self.insert(self.head.prev.key, n)

    def remove(self, key):
        if (self.head.isNil):
            raise EmptyError("ERROR: The route is empty.")
        n = self.query(key)
        if (n.isNil):
            raise KeyNotExistError("ERROR: %s does not exist." % key)
        self.dist += self.tau[n.prev.key, n.next.key] - self.tau[n.prev.key, n.key] - self.tau[n.key, n.next.key]
        if (self.asymFlag):
            self._revDist += self.tau[n.next.key, n.prev.key] - self.tau[n.key, n.prev.key] - self.tau[n.next.key, n.key]
        n.prev.next = n.next
        n.next.prev = n.prev
        if (self.head.key == key):
            self.head = n.next
        self._count -= 1

    def swapNext(self, key):
        n = self.query(key)
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

    def swap(self, keyI, keyJ):
        if (keyI == keyJ):
            raise KeyExistError("ERROR: Cannot swap itself")
        # Old: = = = i j k l m n = = = | -> swap(j, m)
        # New: = = = i m k l j n = = =
        nI = self.query(keyI)
        nJ = self.query(keyJ)
        if (nI.next.key == keyJ):
            return self.swapNext(keyI)
        if (nI.next.next.key == keyJ):
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
        self.insert(self.head.key, n)
        sofarCheapestCost = self.dist if not self.asymFlag else min(self.dist, self._revDist)
        sofarCheapestKey = self.head.key
        if (self._count <= 2):
            return
        cur = self.head
        trvFlag = True
        while (trvFlag):
            self.swapNext(n.key)
            cur = cur.next
            newCost = self.dist if not self.asymFlag else min(self.dist, self._revDist)
            if (newCost < sofarCheapestCost):
                sofarCheapestCost = newCost
                sofarCheapestKey = cur.key
            if (cur.key == self.head.key):
                trvFlag = False
        self.remove(n.key)
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
        return bestRemovalKey

    def impvRemovalReinsert(self):
        # NOTE: This function is useless
        testedFlag = {}
        cur = self.head
        impvedFlag = False
        canImpvFlag = True
        while (canImpvFlag):
            if (cur.key not in testedFlag):
                testedFlag[cur.key] = False            
            if (testedFlag[cur.key] == False):
                canLocImpvFlag = False
                nextNode = cur.next            
                curKey = cur.key
                curPrev = cur.prev.key
                curNext = cur.next.key
                curCopy = cur.clone()
                self.remove(curKey)         
                self.cheapestInsert(curCopy)
                new = self.query(cur.key)
                newPrev = new.prev.key
                newNext = new.next.key
                if (curPrev != newPrev or curNext != newNext):
                    canLocImpvFlag = True
                if (canLocImpvFlag):
                    impvedFlag = True
                    testedFlag[curKey] = False
                    testedFlag[newPrev] = False
                    testedFlag[newNext] = False
                    testedFlag[curPrev] = False
                    testedFlag[curNext] = False
                    cur = nextNode
                else:
                    testedFlag[curKey] = True
                    cur = cur.next
            else:
                cur = cur.next                
            if (min(testedFlag.values()) == True):
                canImpvFlag = False
        return impvedFlag

    def impv2Swap(self):
        # NOTE: This function is useless
        nI = self.head.next
        sofarBestDist = self.dist if not self.asymFlag else min(self.dist, self._revDist)
        while (nI.key != self.head.key):
            # 1. First attempt - 1 swap
            # Old: = = = i j k l m n = = = | -> swap(i, j) 
            # 2. Follow-up attempts
            # New: = = = j i k l m n = = = | -> swap(i, k) -> swap(j, k)
            # New: = = = k j i l m n = = = | -> swap(i, l) -> swap(k, l)
            # New: = = = l j k i m n = = = | -> swap(i, m) -> swap(l, m)
            # New: = = = m j k l i n = = = | -> swap(i, n) -> swap(m, n)
            # 3. Recover to initial status
            # New: = = = n j k l m i = = = | -> swap(i, n)
            # Old: = = = i j k l m n = = = 

            # 1. First attempt
            prevIKey = nI.next.key
            nINextKey = None
            nIKey = nI.key
            nJKey = nI.next.key
            self.swapNext(nI.key)
            newDist = self.dist if not self.asymFlag else min(self.dist, self._revDist)
            if (newDist < sofarBestDist):
                sofarBestDist = newDist
                return True

            # 2. Follow-up attempts
            for _ in range(self._count - 2):
                nINextKey = nI.next.key
                self.swapNext(nI.key)
                self.swap(prevIKey, nINextKey)
                newDist = self.dist if not self.asymFlag else min(self.dist, self._revDist)
                if (newDist < sofarBestDist):
                    sofarBestDist = newDist
                    return True
                prevIKey = nINextKey

            # 3. Recover to initial status
            self.swapNext(nI.key)
            nI = nI.next
        return False

    def impv2Opt(self):
        # NOTE: Again, this function is garbage, I don't understand why the previous-written function runs much faster?
        oriHeadKey = self.head.key
        nI = self.head.next   
        sofarBestDist = self.dist if not self.asymFlag else min(self.dist, self._revDist)
        # print([i.key for i in self.traverse(closeFlag=True)], sofarBestDist, sofarBestDist, "----------------------")
        while (nI.key != self.head.key):
            # 1. First attempt
            # Old: = = = nI nINext nJ nJNext nJ2Next nJ3Next nJ4Next = = = | -> swap(nINext, nJ)        | -> swap(nINext, nJ)
            # 2. Follow-up attempts
            # New: = = = nI nJ nINext nJNext nJ2Next nJ3Next nJ4Next = = = | -> swap(nINext, nJNext)    | -> remove(nJNext)
            #                                                                -> swap(nJ, nJNext)          -> insert(nI, nJNext)
            # New: = = = nI nJNext nJ nINext nJ2Next nJ3Next nJ4Next = = = | -> swap(nINext, n2JNext)   | -> remove(nJ2Next)
            #                                                                -> swap(nJ, nJ2Next)         -> insert(nI, nJ2Next)
            #                                                                -> swap(nJNext, nJ2Next)
            # New: = = = nI nJ2Next nJNext nJ nINext nJ3Next nJ4Next = = = | -> swap(nINext, n3JNext)   | -> remove(nJ3Next)
            #                                                                -> swap(nJ, nJ3Next)         -> insert(nI, nJ3Next)
            #                                                                -> swap(nJNext, nJ3Next) 
            #                                                                -> swap(nJ2Next, nJ3Next)
            # New: = = = nI nJ3Next nJ2Next nJNext nJ nINext nJ4Next = = =
            # ...
            # 3. Recover to initial status (Until: nJXNext == nIPrev)
            # Old: nI nJ = = = nINext nJNext (nIPrev) | -> swap(nJNext, nI)
            # New: nI nJNext nJ = = = nINext
            # self.reverse()

            # 1. First attempt
            nJ = nI.next.next
            nINext = nI.next
            self.swapNext(nINext.key)
            newDist = self.dist if not self.asymFlag else min(self.dist, self._revDist)
            # print([i.key for i in self.traverse(closeFlag=True)], newDist, sofarBestDist, "nI", nI.key, "SwapNext: ", nINext.next.key)
            if (newDist < sofarBestDist):
                sofarBestDist = newDist
                self.rehead(oriHeadKey)
                return True

            # 2. Follow-up attempts
            for _ in range(self._count - 4):
                nJXNext = nINext.next
                nJXNextClone = nJXNext.clone()
                self.remove(nJXNext.key)
                self.insert(nI.key, nJXNextClone)
                newDist = self.dist if not self.asymFlag else min(self.dist, self._revDist)
                # print([i.key for i in self.traverse(closeFlag=True)], newDist, sofarBestDist, "nI", nI.key, "RemoveAndInsert: ", nJXNext.key)
                if (newDist < sofarBestDist):
                    sofarBestDist = newDist
                    self.rehead(oriHeadKey)
                    return True

            # 3. Recover to intial status
            nJXNext = nINext.next
            self.swapNext(nJXNext.key)
            self.reverse()
            self.rehead(oriHeadKey)
            # print([i.key for i in self.traverse(closeFlag=True)], newDist, sofarBestDist, "nI", nI.key, "Recover")
            nI = nI.next
        return False

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
        return self.root == self.nil

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
