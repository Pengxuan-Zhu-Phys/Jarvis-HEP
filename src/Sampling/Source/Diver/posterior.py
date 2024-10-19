#!/usr/bin/env python3 

import numpy as np

class Point:
    def __init__(self, vector=None, weight=None):
        self.vector = vector
        self.weight = weight
        self.next = None

class Node:
    def __init__(self, lowerbounds=None, upperbounds=None):
        self.trunk = None
        self.branchA = None
        self.branchB = None
        self.firstnewpt = None
        self.lastnewpt = None
        self.population = 0.0
        self.newpopulation = 0
        self.branchesDifferInDim = 0
        self.lowerbounds = np.array(lowerbounds) if lowerbounds is not None else None
        self.upperbounds = np.array(upperbounds) if upperbounds is not None else None

class ListNode:
    def __init__(self):
        self.thisnode = None
        self.next = None
        self.prev = None

class Duplicate:
    def __init__(self):
        self.next = None
        self.point = None
        self.stoppingNode = None

class PostModule:
    def __init__(self):
        self.root = None
        self.firstListNode = None
        self.currentListNode = None
        self.firstDuplicate = None
        self.lastDuplicate = None
        self.ranges = None
        self.D = 0
        self.totalCells = 0
        self.debug = False
        self.maxNodePop = None
        self.duplicate_tol = 10.0

    def ini_tree(self, lowerbounds, upperbounds, maxpop):
        """Initialize the root of the tree and the list of tree nodes."""
        self.totalCells = 1
        self.D = len(lowerbounds)
        self.root = Node(lowerbounds, upperbounds)
        self.ranges = upperbounds - lowerbounds
        self.firstListNode = ListNode()
        self.maxNodePop = maxpop

    def clear_tree(self):
        """Clear the tree."""
        self.totalCells = 1
        self.burn_tree(self.root)

    def burn_tree(self, currentNode):
        """Recursively burns off the tree from the currentNode."""
        if currentNode.branchesDifferInDim != 0:
            self.burn_tree(currentNode.branchA)
            self.burn_tree(currentNode.branchB)
            currentNode.branchA = None
            currentNode.branchB = None

    def get_weight(self, vector, prior, context):
        """Find the posterior weight of a single individual using the tree."""
        individual = Point(vector=vector, weight=0)
        self.climb_tree(individual, self.root, justLooking=True)
        return individual.weight * self.totalCells * prior(vector, context)

    def grow_tree(self, X, prior, context):
        """Grow the tree with points in a new generation."""
        NP = len(X['values'])
        self.currentListNode = self.firstListNode
        self.currentListNode.next = None

        if self.debug:
            print('Current number of cells:', self.totalCells)
            print('About to climb')

        for i in range(NP):
            individual = Point(vector=X['vectors'][i], weight=X['weights'][i])
            self.climb_tree(individual, self.root)

        self.tend_tree(self.firstListNode)

        while self.firstDuplicate is not None:
            if self.debug:
                print('Dealing with duplicates...')
            self.currentListNode = self.firstListNode
            self.currentListNode.next = None
            self.straggle_up_tree(self.firstDuplicate)
            self.tend_tree(self.firstListNode)

        X['weights'] *= self.totalCells * np.array([prior(vec, context) for vec in X['vectors']])

    def tend_tree(self, workingListNode):
        """Check and grow the tree at each node in the linked list."""
        while workingListNode is not None:
            if workingListNode.prev is None:
                workingListNode = workingListNode.next
                continue

            node = workingListNode.thisnode

            if node.population > self.maxNodePop:
                self.grow_branches(node)
            else:
                nodeWeight = np.prod(node.upperbounds - node.lowerbounds)
                individual = node.firstnewpt
                for _ in range(node.newpopulation):
                    individual.weight = nodeWeight
                    individual = individual.next

            nextNode = workingListNode.next
            workingListNode = None  # Deallocate this node
            workingListNode = nextNode

    def straggle_up_tree(self, straggler):
        """Send stragglers climbing up the tree after initial points."""
        self.climb_tree(straggler.point, straggler.stoppingNode)
        if straggler.next is not None:
            self.straggle_up_tree(straggler.next)
        straggler = None  # Deallocate straggler

    def climb_tree(self, individual, currentNode, justLooking=False):
        """Climb a level in the binary space partition tree."""
        if currentNode.branchesDifferInDim == 0:
            if justLooking:
                individual.weight = np.prod(currentNode.upperbounds - currentNode.lowerbounds)
            else:
                currentNode.population += 1
                self.add_to_end_of_pt_list(individual, currentNode)

                if currentNode.firstnewpt.next is None:
                    newNode = ListNode()
                    newNode.prev = self.currentListNode
                    self.currentListNode.next = newNode
                    self.currentListNode = newNode
                    self.currentListNode.thisnode = currentNode
        else:
            if individual.vector[currentNode.branchesDifferInDim] < currentNode.branchA.upperbounds[currentNode.branchesDifferInDim]:
                self.climb_tree(individual, currentNode.branchA, justLooking)
            else:
                self.climb_tree(individual, currentNode.branchB, justLooking)

    def add_to_end_of_pt_list(self, individual, currentNode):
        """Add a point to the end of the linked list of points."""
        if currentNode.newpopulation == 0:
            currentNode.firstnewpt = individual
            currentNode.lastnewpt = individual
        else:
            currentNode.lastnewpt.next = individual
            currentNode.lastnewpt = individual
        currentNode.newpopulation += 1

    def grow_branches(self, currentNode):
        """Add branches to a node, and sort the population into them."""
        A = Node(currentNode.lowerbounds, currentNode.upperbounds)
        B = Node(currentNode.lowerbounds, currentNode.upperbounds)
        currentNode.branchA = A
        currentNode.branchB = B
        A.upperbounds[currentNode.branchesDifferInDim] = (currentNode.upperbounds[currentNode.branchesDifferInDim] + currentNode.lowerbounds[currentNode.branchesDifferInDim]) / 2
        B.lowerbounds[currentNode.branchesDifferInDim] = A.upperbounds[currentNode.branchesDifferInDim]
        currentNode.population = 0
        currentNode.newpopulation = 0
        self.totalCells += 1
