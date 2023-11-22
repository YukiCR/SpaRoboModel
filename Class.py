"""
providing the class needed
"""
import numpy as np
import Transform as mytf
from typing import List, Dict, Tuple

class PoseNode:
    """
    ### PoseNode
    a Tree class for dealing with pose\n
    add and remove parent&child by addChild and rmChild method
    """
    def __init__(self, name:str="UnNamed", parent=None):
        """
        do some init \n
        name: string \n
        parent: an instance of PoseNode \n
        name, parent, children are properties\
            that manages __name, __parent and\
            __children 
        """
        self.name = name
        self.parent = None
        self.__children = []

        if isinstance(parent, PoseNode):
            parent.addChild(self)

    @property
    def name(self):
        return self.__name
    
    @name.setter
    def name(self, name):
        if isinstance(name, str):
            self.__name = name
        else:
            raise TypeError("name should be a string")
        

    @property
    def parent(self):
        """
        parent property \n
        getter: returns self.__parent \n
        setter: make sure type(parent) can only be PoseNode \
            or None
        """
        return self.__parent
    
    @parent.setter
    def parent(self, parent):
        if isinstance(parent, PoseNode) or parent == None:
            self.__parent = parent
        else:
            raise TypeError("parent should be an instance of PoseNode or None")
        # raise ValueError("parent can only be modified by add and rm")
        

    @property
    def children(self):
        return self.__children
    
    @children.setter
    def children(self, children):
        if isinstance(children, PoseNode):
            self.__children = children
            print("be careful when using children setter,\
                  recommand addChild and rmChild")
        else:
            raise TypeError("parent should be an instance of PoseNode")
        # raise ValueError("children can only modified by add and rm")
    

    def addChild(self, child):
        """
        add a child for self, it would automatically add parent for the child as well
        """
        if isinstance(child, PoseNode):
            if child not in self.__children: # avoid repetitive append
                self.__children.append(child)
            child.parent = self
        else:
            raise TypeError("child should be an instance of PoseNode")


    def rmChild(self, child):
        if isinstance(child, PoseNode):
            self.__children.remove(child)
            child.parent = None
        else:
            raise TypeError("child should be an instance of PoseNode")
        
    def print(self, level=0):
        """
        a print written by GPT, print recursively
        """
        if self.parent != None:
            print("\t" * level + self.name + "(" + self.parent.name + ")")
        else:
            print("\t" * level + self.name)
        for child in self.children:
            child.print_tree(level + 1)





class joint(PoseNode):
    """
    ### Joint
    a joint class storing the mass and pose information \n
    it derives the PoseNode class, so parent and children can be appointed\n
    + m:mass, a,b and I are a_ii, b_ii and I_ii\n
        p is position p(i,i-1), R is quaternion R(i,i-1)\n
        using Quaternion takes less storage, but R computes faster\n
        dp,dtheta: diravative of p and theta, w: angle velocity
    + both base and joint can use this class, as for base, theta is a constant(0)\n
        as for joint, p and q are constant
    + T = [R,P;0,1], T = T(p)@T(R)@T(theta)
    """
    def __init__(self, name:str="UnNamed", parent=None, m:float=0, \
                 a:np.ndarray=np.zeros(3), b:np.ndarray=np.zeros(3), I:np.ndarray=np.zeros((3,3)), \
                 p:np.ndarray=np.zeros(3), R:np.ndarray=np.eye(3), theta:float=0,\
                 dp:np.ndarray=np.zeros(3), w:np.ndarray=np.zeros(3), dtheta:float=0):
        super().__init__(name, parent)
        self.__m = m
        self.__a = a
        self.__b = b
        self.__I = I
        self.p = p
        self.R = R
        self.theta = theta
        self.dp = dp
        self.w = w
        self.dtheta = dtheta

    @property
    def T(self):
        """
        The T matrix form i-1 to i \n
        every time call joint.T, the getter would actively calculate the T mat
        """
        # T(p)@T(R)@T(theta)
        return  mytf.augment(R=self.R, p=self.p) @ mytf.augment(R=mytf.Rz(self.theta))
    
    @T.setter
    def T(self,input):
        raise ValueError("T is read only")
    
    
    




class robot:
    def __init__(self) -> None:
        """
        by putting joints into robots, more info can be indicated \n
        r: position of center fo mass of joint \n
        T: the T mat T(i,0), the T property of joint is T(i,i-1) \n
        k: Z_hat(i,0) \n
        """
        self.__joints:List[joint] = [] # joints series, do not support parrial robot 
        self.__r:List[np.ndarray] = {} # position of center fo mass of joint
        self.__T:List[np.ndarray] = {} # T(i,0)
        self.__k:List[np.ndarray] = {} # Z_hat(i,0) 
    
    def addJoint(self, child:joint, parent:joint):
        if isinstance(child, joint) and isinstance(parent, joint):
            parent.addChild(child) # link child to parent
        self.__joints.append(child) # append child to the joint list
        self.update()

    def update(self):
        # to be fixed
        for i,joint in enumerate(self.__joints):
            if joint.parent == None:
                self.__T[i] = joint.T
            else:
                self.__T[i] = self.__T[i-1] @ joint.T
            self.__k[i] = self.__T[i][0:3,-1]
            self.__r[i] = self.__T[i] @ ......
            
            
                

    






    
    
    






    
# check cross
r = robot()
a = joint()
b = joint()
r.addJoint(b,a)

    