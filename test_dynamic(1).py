import numpy as np
import matplotlib.pyplot as plt

def system_dynamics(state, t,u):

    # 本体位姿，  6个theta， 本体位姿一阶导，   Dtheta,     
    x0, y0, z0, q0, q1, q2, q3, theta1, theta2, theta3, theta4, theta5, theta6, vx0, vy0, vz0, wx0, wy0, wz0, dtheta1, dtheta2, dtheta3, dtheta4, dtheta5, dtheta6 = state
#动力学所需参数：theta1-6,q0-3,x0-z0
# dtheta1-6,vx0-vz0,wx0-wz0

    w0 = np.array([[wx0],[wy0],[wz0]])
    q = np.array([[q0],[q1],[q2],[q3]])

    """
    T = [R p]
        [0 1]
    """
    def rot_x(x):
        Tx = np.array([[1,0,0,0],[0,np.cos(x),-np.sin(x),0],[0,np.sin(x),np.cos(x),0],[0,0,0,1]])
        return Tx
    def rot_y(y):
        Ty = np.array([[np.cos(y),0,np.sin(y),0],[0,1,0,0],[-np.sin(y),0,np.cos(y),0],[0,0,0,1]])
        return Ty
    def rot_z(z):
        Tz = np.array([[np.cos(z),-np.sin(z),0,0],[np.sin(z),np.cos(z),0,0],[0,0,1,0],[0,0,0,1]])
        return Tz
    
    def tran_x(x):
        Tx = np.array([[1,0,0,x],[0,1,0,0],[0,0,1,0],[0,0,0,1]])
        return Tx
    def tran_y(y):
        Ty = np.array([[1,0,0,0],[0,1,0,y],[0,0,1,0],[0,0,0,1]])
        return Ty
    def tran_z(z):
        Tz = np.array([[1,0,0,0],[0,1,0,0],[0,0,1,z],[0,0,0,1]])
        return Tz
    
    def ccsz(a):
        """cross"""
        T = np.array([[0,-a[2][0],a[1][0]],[a[2][0],0,-a[0][0]],[-a[1][0],a[0][0],0]])
        return T
    def rot_sys(q):
        """Q to T"""
        H = np.array([[1-2*q[2][0]**2-2*q[3][0]**2 , 2*q[1][0]*q[2][0]-2*q[0][0]*q[3][0] , 2*q[1][0]*q[3][0]+2*q[0][0]*q[2][0] , 0], 
                    [2*q[1][0]*q[2][0]+2*q[0][0]*q[3][0] , 1-2*q[1][0]**2-2*q[3][0]**2 , 2*q[2][0]*q[3][0]-2*q[0][0]*q[1][0] , 0],
                    [2*q[1][0]*q[3][0]-2*q[0][0]*q[2][0] , 2*q[2][0]*q[3][0]+2*q[0][0]*q[1][0] , 1-2*q[1][0]**2-2*q[2][0]**2 , 0],
                    [0 , 0 , 0 , 1]])
        return H

    #各坐标系变换
    # I T
    # 0
    # [P' 1]'
    TI0 = tran_x(x0) @ tran_y(y0) @ tran_z(z0) @ rot_sys(q)
    T01 = tran_x(0.3570) @ tran_y(-0.0095) @ tran_z(0.419) @ rot_z(-np.pi/2) @ rot_z(theta1)
    T12 = tran_z(0.3) @ rot_x(np.pi/2) @ rot_y(np.pi/2) @ rot_z(theta2)
    T23 = tran_z(-0.3) @ tran_x(0.83) @ rot_z(np.pi) @ rot_z(theta3)
    T34 = tran_z(0.3) @ rot_z(np.pi/2) @ rot_x(np.pi/2) @ rot_z(theta4)
    T45 = tran_z(-0.7) @ rot_x(np.pi/2) @ rot_y(np.pi) @ rot_z(theta5)
    T56 = tran_y(-0.1) @ rot_x(np.pi/2) @ rot_z(theta6)
    T6E = tran_z(0.2345)
    
    # T, P=0
    T_I0 = rot_sys(q)
    T_01 = rot_z(-np.pi/2) @ rot_z(theta1)
    T_12 = rot_x(np.pi/2) @ rot_y(np.pi/2) @ rot_z(theta2)
    T_23 = rot_z(np.pi) @ rot_z(theta3)
    T_34 = rot_z(np.pi/2) @ rot_x(np.pi/2) @ rot_z(theta4)
    T_45 = rot_x(np.pi/2) @ rot_y(np.pi) @ rot_z(theta5)
    T_56 = rot_x(np.pi/2) @ rot_z(theta6)

    T_I03 = T_I0[0:3,0:3]
    T_013 = T_01[0:3,0:3]
    T_123 = T_12[0:3,0:3]
    T_233 = T_23[0:3,0:3]
    T_343 = T_34[0:3,0:3]
    T_453 = T_45[0:3,0:3]
    T_563 = T_56[0:3,0:3]

    # 第k个连体坐标系的Z轴的方向向量在惯性系下的描述
    k1 = T_I0 @ T_01 @ np.array([[0],[0],[1],[1]])
    k2 = T_I0 @ T_01 @ T_12 @ np.array([[0],[0],[1],[1]])
    k3 = T_I0 @ T_01 @ T_12 @ T_23 @ np.array([[0],[0],[1],[1]])
    k4 = T_I0 @ T_01 @ T_12 @ T_23 @ T_34 @ np.array([[0],[0],[1],[1]])
    k5 = T_I0 @ T_01 @ T_12 @ T_23 @ T_34 @ T_45 @ np.array([[0],[0],[1],[1]])
    k6 = T_I0 @ T_01 @ T_12 @ T_23 @ T_34 @ T_45 @ T_56 @ np.array([[0],[0],[1],[1]])

    k1 = k1[0:3,:]
    k2 = k2[0:3,:]
    k3 = k3[0:3,:]
    k4 = k4[0:3,:]
    k5 = k5[0:3,:]
    k6 = k6[0:3,:]

    # 质心位置在惯性系下表述
    r0 = TI0 @ np.array([[0],[0],[0],[1]])
    r1 = TI0 @ T01 @ np.array([[0],[0],[0.15],[1]])
    r2 = TI0 @ T01 @ T12 @ np.array([[0.2702],[0],[-0.2513],[1]])
    r3 = TI0 @ T01 @ T12 @ T23 @ np.array([[0],[0],[0.15],[1]])
    r4 = TI0 @ T01 @ T12 @ T23 @ T34 @ np.array([[0],[0],[-0.35],[1]])
    r5 = TI0 @ T01 @ T12 @ T23 @ T34 @ T45 @ np.array([[0],[-0.0338],[0],[1]])
    r6 = TI0 @ T01 @ T12 @ T23 @ T34 @ T45 @ T56 @ np.array([[0],[0],[0.0750],[1]])

    # 连体坐标系原点在惯性系下的表述
    p1 = TI0 @ T01 @ np.array([[0],[0],[0],[1]])
    p2 = TI0 @ T01 @ T12 @ np.array([[0],[0],[0],[1]])
    p3 = TI0 @ T01 @ T12 @ T23 @ np.array([[0],[0],[0],[1]])
    p4 = TI0 @ T01 @ T12 @ T23 @ T34 @ np.array([[0],[0],[0],[1]])
    p5 = TI0 @ T01 @ T12 @ T23 @ T34 @ T45 @ np.array([[0],[0],[0],[1]])
    p6 = TI0 @ T01 @ T12 @ T23 @ T34 @ T45 @ T56 @ np.array([[0],[0],[0],[1]])
    pe = TI0 @ T01 @ T12 @ T23 @ T34 @ T45 @ T56 @ T6E @ np.array([[0],[0],[0],[1]])

    r0 = r0[0:3,:]
    r1 = r1[0:3,:]
    r2 = r2[0:3,:]
    r3 = r3[0:3,:]
    r4 = r4[0:3,:]
    r5 = r5[0:3,:]
    r6 = r6[0:3,:]
    p1 = p1[0:3,:]
    p2 = p2[0:3,:]
    p3 = p3[0:3,:]
    p4 = p4[0:3,:]
    p5 = p5[0:3,:]
    p6 = p6[0:3,:]
    pe = pe[0:3,:]

    m0 = 400
    m1 = 6
    m2 = 5
    m3 = 5
    m4 = 4
    m5 = 3
    m6 = 2

    # 连体系下的惯量， 原点是？
    I00 = np.array([[30,0.26,0.37],[0.26,28,-0.29],[0.37,-0.29,32]])
    I11 = np.array([[0.15,0,0],[0,0.15,0],[0,0,0.075]])
    I22 = np.array([[0.0926,0,0.1315],[0,0.9053,0],[0.1315,0,0.8451]])
    I33 = np.array([[0.105,0,0],[0,0.105,0],[0,0,0.0294]])
    I44 = np.array([[0.2498,0,0],[0,0.2498,0],[0,0,0.0196]])
    I55 = np.array([[0.0330,0,0],[0,0.0172,0],[0,0,0.0260]])
    I66 = np.array([[0.05152,0,0],[0,0.05152,0],[0,0,0.02192]])

    # 惯性系下的惯量，
    I0I = T_I03 @ I00 @ np.transpose(T_I03)
    I1I = T_I03 @ T_013 @ I11 @ np.transpose(T_I03 @ T_013)
    I2I = T_I03 @ T_013 @ T_123 @ I22 @ np.transpose(T_I03 @ T_013 @ T_123)
    I3I = T_I03 @ T_013 @ T_123 @ T_233 @ I33 @ np.transpose(T_I03 @ T_013 @ T_123 @ T_233)
    I4I = T_I03 @ T_013 @ T_123 @ T_233 @ T_343 @ I44 @ np.transpose(T_I03 @ T_013 @ T_123 @ T_233 @ T_343)
    I5I = T_I03 @ T_013 @ T_123 @ T_233 @ T_343 @ T_453 @ I55 @ np.transpose(T_I03 @ T_013 @ T_123 @ T_233 @ T_343 @ T_453)
    I6I = T_I03 @ T_013 @ T_123 @ T_233 @ T_343 @ T_453 @ T_563 @ I66 @ np.transpose(T_I03 @ T_013 @ T_123 @ T_233 @ T_343 @ T_453 @ T_563)

    r01 = r1 - r0
    r02 = r2 - r0
    r03 = r3 - r0
    r04 = r4 - r0
    r05 = r5 - r0
    r06 = r6 - r0
    M = m0 + m1 + m2 + m3 + m4 + m5 + m6
    rs =(r0*m0 + r1*m1 + r2*m2 + r3*m3 + r4*m4 + r5*m5 + r6*m6)/M

    #算JTI
    JT_0 = np.array([[0],[0],[0]])
    JT_11 = ccsz(k1)@(r1-p1)
    JT_21 = ccsz(k1)@(r2-p1)
    JT_22 = ccsz(k2)@(r2-p2)
    JT_31 = ccsz(k1)@(r3-p1)
    JT_32 = ccsz(k2)@(r3-p2)
    JT_33 = ccsz(k3)@(r3-p3)
    JT_41 = ccsz(k1)@(r4-p1)
    JT_42 = ccsz(k2)@(r4-p2)
    JT_43 = ccsz(k3)@(r4-p3)
    JT_44 = ccsz(k4)@(r4-p4)
    JT_51 = ccsz(k1)@(r5-p1)
    JT_52 = ccsz(k2)@(r5-p2)
    JT_53 = ccsz(k3)@(r5-p3)
    JT_54 = ccsz(k4)@(r5-p4)
    JT_55 = ccsz(k5)@(r5-p5)
    JT_61 = ccsz(k1)@(r6-p1)
    JT_62 = ccsz(k2)@(r6-p2)
    JT_63 = ccsz(k3)@(r6-p3)
    JT_64 = ccsz(k4)@(r6-p4)
    JT_65 = ccsz(k5)@(r6-p5)
    JT_66 = ccsz(k6)@(r6-p6)
    JT1 = np.concatenate((JT_11,JT_0,JT_0,JT_0,JT_0,JT_0),axis=1)
    JT2 = np.concatenate((JT_21,JT_22,JT_0,JT_0,JT_0,JT_0),axis=1)
    JT3 = np.concatenate((JT_31,JT_32,JT_33,JT_0,JT_0,JT_0),axis=1)
    JT4 = np.concatenate((JT_41,JT_42,JT_43,JT_44,JT_0,JT_0),axis=1)
    JT5 = np.concatenate((JT_51,JT_52,JT_53,JT_54,JT_55,JT_0),axis=1)
    JT6 = np.concatenate((JT_61,JT_62,JT_63,JT_64,JT_65,JT_66),axis=1)

    #算JR
    JR_1 = k1
    JR_2 = k2
    JR_3 = k3
    JR_4 = k4
    JR_5 = k5
    JR_6 = k6
    JR1 = np.concatenate((JR_1,JT_0,JT_0,JT_0,JT_0,JT_0),axis=1)
    JR2 = np.concatenate((JR_1,JR_2,JT_0,JT_0,JT_0,JT_0),axis=1)
    JR3 = np.concatenate((JR_1,JR_2,JR_3,JT_0,JT_0,JT_0),axis=1)
    JR4 = np.concatenate((JR_1,JR_2,JR_3,JR_4,JT_0,JT_0),axis=1)
    JR5 = np.concatenate((JR_1,JR_2,JR_3,JR_4,JR_5,JT_0),axis=1)
    JR6 = np.concatenate((JR_1,JR_2,JR_3,JR_4,JR_5,JR_6),axis=1)


    mat11 = np.array([[M,0,0],[0,M,0],[0,0,M]])
    mat12 = M * np.transpose(ccsz(rs-r0))
    mat13 = m1*JT1 + m2*JT2 + m3*JT3 + m4*JT4 + m5*JT5 + m6*JT6
    mat21 = np.transpose(mat12)
    mat22 = I0I + I1I + I2I + I3I + I4I + I5I + I6I - m1 *ccsz(r01)@ccsz(r01) - m2 *ccsz(r02)@ccsz(r02) - m3 *ccsz(r03)@ccsz(r03) - m4 *ccsz(r04)@ccsz(r04) - m5 *ccsz(r05)@ccsz(r05) - m6 *ccsz(r06)@ccsz(r06)
    mat23 = I1I @ JR1 + I2I @ JR2 + I3I @ JR3 + I4I @ JR4 + I5I @ JR5 + I6I @ JR6 + m1*ccsz(r01)@JT1 + m2*ccsz(r02)@JT2 + m3*ccsz(r03)@JT3 + m4*ccsz(r04)@JT4 + m5*ccsz(r05)@JT5 + m6*ccsz(r06)@JT6
    mat31 = np.transpose(mat13)
    mat32 = np.transpose(mat23)
    mat33 = np.transpose(JR1)@I1I@JR1 + np.transpose(JR2)@I2I@JR2 + np.transpose(JR3)@I3I@JR3 + np.transpose(JR4)@I4I@JR4 + np.transpose(JR5)@I5I@JR5 + np.transpose(JR6)@I6I@JR6 + m1*np.transpose(JT1)@JT1 + m2*np.transpose(JT2)@JT2 + m3*np.transpose(JT3)@JT3 + m4*np.transpose(JT4)@JT4 + m5*np.transpose(JT5)@JT5 + m6*np.transpose(JT6)@JT6
    H1 = np.concatenate((mat11,mat12,mat13),axis=1)
    H2 = np.concatenate((mat21,mat22,mat23),axis=1)
    H3 = np.concatenate((mat31,mat32,mat33),axis=1)
    H = np.concatenate((H1,H2,H3),axis=0)
    







    #计算非线性项， C

    R10 = (rot_z(-theta1)@rot_z(np.pi/2))[0:3,0:3]
    R21 = (rot_z(-theta2)@rot_y(-np.pi/2) @ rot_x(-np.pi/2))[0:3,0:3]
    R32 = (rot_z(-theta3)@rot_z(np.pi))[0:3,0:3]
    R43 = (rot_z(-theta4)@rot_z(-np.pi/2)@rot_y(-np.pi/2))[0:3,0:3]
    R54 = (rot_z(-theta5)@rot_z(np.pi)@rot_x(np.pi/2))[0:3,0:3]
    R65 = (rot_z(-theta6)@rot_x(-np.pi/2))[0:3,0:3]

    R01 = (rot_z(-np.pi/2)@rot_z(theta1))[0:3,0:3]
    R12 = (rot_x(np.pi/2) @ rot_y(np.pi/2)@rot_z(theta2))[0:3,0:3]
    R23 = (rot_z(np.pi)@rot_z(theta3))[0:3,0:3]
    R34 = (rot_y(np.pi/2)@rot_z(np.pi/2)@rot_z(theta4))[0:3,0:3]
    R45 = (rot_x(-np.pi/2)@rot_z(np.pi)@rot_z(theta5))[0:3,0:3]
    R56 = (rot_x(np.pi/2)@rot_z(theta6))[0:3,0:3]

    P01 = np.array([[0.3570],[-0.0095],[0.419]])
    Pc11 = np.array([[0],[0],[0.15]])
    P12 = np.array([[0],[0],[0.3]])
    Pc22 = np.array([[0.2702],[0],[-0.2513]])
    P23 = np.array([[0.83],[0],[-0.3]])
    Pc33 = np.array([[0],[0],[0.15]])
    P34 = np.array([[0],[0],[0.3]])
    Pc44 = np.array([[0],[0],[-0.35]])
    P45 = np.array([[0],[0],[-0.7]])
    Pc55 = np.array([[0],[-0.0338],[0]])
    P56 = np.array([[0],[-0.1],[0]])
    Pc66 = np.array([[0],[0],[0.075]])

    I1C = I11
    I2C = I22
    I3C = I33
    I4C = I44
    I5C = I55
    I6C = I66

    w11 = R10 @ w0 + dtheta1 * np.array([[0],[0],[1]])
    dw11 = R10 @ccsz(w0) @ (dtheta1 * np.array([[0],[0],[1]]))
    dv11 = R10 @ (ccsz(w0) @ (ccsz(w0)@P01))
    dvc11 = ccsz(dw11) @ Pc11 + ccsz(w11) @ (ccsz(w11) @ Pc11) + dv11
    F11 = m1 * dvc11
    N11 = I1C @ dw11 + ccsz(w11) @ I1C @ w11

    w22 = R21 @ w11 + dtheta2 * np.array([[0],[0],[1]])
    dw22 = R21 @ dw11 + R21 @ ccsz(w11) @ (dtheta2 * np.array([[0],[0],[1]]))
    dv22 = R21 @ (ccsz(dw11)@P12 + ccsz(w11) @ (ccsz(w11) @ P12) + dv11)
    dvc22 = ccsz(dw22) @ Pc22 + ccsz(w22) @ (ccsz(w22) @ Pc22) + dv22
    F22 = m2 * dvc22
    N22 = I2C @ dw22 + ccsz(w22) @ I2C @ w22

    w33 = R32 @ w22 + dtheta3 * np.array([[0],[0],[1]])
    dw33 = R32 @ dw22 + R32 @ ccsz(w22) @ (dtheta3 * np.array([[0],[0],[1]]))
    dv33 = R32 @ (ccsz(dw22)@P23 + ccsz(w22) @ (ccsz(w22) @ P23) + dv22)
    dvc33 = ccsz(dw33) @ Pc33 + ccsz(w33) @ (ccsz(w33) @ Pc33) + dv33
    F33 = m3 * dvc33
    N33 = I3C @ dw33 + ccsz(w33) @ I3C @ w33

    w44 = R43 @ w33 + dtheta4 * np.array([[0],[0],[1]])
    dw44 = R43 @ dw33 + R43 @ ccsz(w33) @ (dtheta4 * np.array([[0],[0],[1]]))
    dv44 = R43 @ ((ccsz(dw33)@P34) + ccsz(w33) @ (ccsz(w33) @ P34) + dv33)
    dvc44 = ccsz(dw44) @ Pc44 + ccsz(w44) @ (ccsz(w44) @ Pc44) + dv44
    F44 = m4 * dvc44
    N44 = I4C @ dw44 + ccsz(w44) @ I4C @ w44

    w55 = R54 @ w44 + dtheta5 * np.array([[0],[0],[1]])
    dw55 = R54 @ dw44 + R54 @ ccsz(w44) @ (dtheta5 * np.array([[0],[0],[1]]))
    dv55 = R54 @ ((ccsz(dw44)@P45) + ccsz(w44) @ (ccsz(w44) @ P45) + dv44)
    dvc55 = ccsz(dw55) @ Pc55 + ccsz(w55) @ (ccsz(w55) @ Pc55) + dv55
    F55 = m5 * dvc55
    N55 = I5C @ dw55 + ccsz(w55) @ I5C @ w55

    w66 = R65 @ w55 + dtheta6 * np.array([[0],[0],[1]])
    dw66 = R65 @ dw55 + R65 @ ccsz(w55) @ (dtheta6 * np.array([[0],[0],[1]]))
    dv66 = R65 @ ((ccsz(dw55)@P56) + ccsz(w66) @ (ccsz(w55) @ P56) + dv55)
    dvc66 = ccsz(dw66) @ Pc66 + ccsz(w66) @ (ccsz(w66) @ Pc66) + dv66
    F66 = m6 * dvc66
    N66 = I6C @ dw66 + ccsz(w66) @ I6C @ w66

    f66 = F66
    n66 = N66 + ccsz(Pc66) @ F66 
    c6 = np.transpose(n66) @ np.array([[0],[0],[1]])

    f55 = F55 + R56@f66
    n55 = N55 + R56@n66 + ccsz(Pc55) @ F55 + ccsz(P56) @ R56 @ f66
    c5 = np.transpose(n55) @ np.array([[0],[0],[1]])

    f44 = F44 + R45@f55 
    n44 = N44 + R45@n55 + ccsz(Pc44) @ F44 + ccsz(P45) @ R45 @ f55
    c4 = np.transpose(n44) @ np.array([[0],[0],[1]])

    f33 = F33 + R34@f44 
    n33 = N33 + R34@n44 + ccsz(Pc33) @ F33 + ccsz(P34) @ R34 @ f44
    c3 = np.transpose(n33) @ np.array([[0],[0],[1]])

    f22 = F22 + R23@f33 
    n22 = N22 + R23@n33 + ccsz(Pc22) @ F22 + ccsz(P23) @ R23 @ f33
    c2 = np.transpose(n22) @ np.array([[0],[0],[1]])

    f11 = F11 + R12@f22 
    n11 = N11 + R12@n22 + ccsz(Pc11) @ F11 + ccsz(P12) @ R12 @ f22
    c1 = np.transpose(n11) @ np.array([[0],[0],[1]])

    f00 = R01@f11 
    n00 = R01@n11 + ccsz(P01) @ R01 @ f11 
    
    c_1 = f00[0][0]
    c_2 = f00[1][0]
    c_3 = f00[2][0]
    c_4 = n00[0][0]
    c_5 = n00[1][0]
    c_6 = n00[2][0]
    C = np.array([[c_1],[c_2],[c_3],[c_4],[c_5],[c_6],[c1[0][0]],[c2[0][0]],[c3[0][0]],[c4[0][0]],[c5[0][0]],[c6[0][0]]])
    print(C)
    d = np.array([[0],[0],[0],[0],[0],[0],[0.05*np.cos(0.05*np.pi*t)],[0.05*np.cos(0.05*np.pi*t)],[0.05*np.cos(0.05*np.pi*t)],[0.005*np.cos(0.05*np.pi*t)],[0.002*np.cos(0.05*np.pi*t)],[0.001*np.cos(0.05*np.pi*t)]])
    ddqdt = np.linalg.inv(H)@(d)

    ddq1 = ddqdt[0][0]
    ddq2 = ddqdt[1][0]
    ddq3 = ddqdt[2][0]
    ddq4 = ddqdt[3][0]
    ddq5 = ddqdt[4][0]
    ddq6 = ddqdt[5][0]
    ddq7 = ddqdt[6][0]
    ddq8 = ddqdt[7][0]
    ddq9 = ddqdt[8][0]
    ddq10 = ddqdt[9][0]
    ddq11 = ddqdt[10][0]
    ddq12 = ddqdt[11][0]
    qv = np.array([[q1],[q2],[q3]])
    dq = 0.5*(np.concatenate((-np.transpose(qv),(np.array([[q0,0,0],[0,q0,0],[0,0,q0]])-ccsz(qv))),axis=0))@w0
    dq0dt = dq[0][0]
    dq1dt = dq[1][0]
    dq2dt = dq[2][0]
    dq3dt = dq[3][0]
    
    return [vx0, vy0, vz0, dq0dt, dq1dt, dq2dt, dq3dt, dtheta1, dtheta2, dtheta3, dtheta4, dtheta5, dtheta6, ddq1, ddq2, ddq3, ddq4, ddq5, ddq6, ddq7, ddq8, ddq9, ddq10, ddq11, ddq12]


def controller(state):
    """
    ---
    """
    control = np.array([[0],[0],[0],[0],[0],[0],[0],[0],[0],[0],[0],[0]])
    return control.ravel()

def runge_kutta_step(func, t, state, u, dt):
    k1 = np.array(func(state, t, u))
    k2 = np.array(func(state + 0.5 * dt * k1, t + 0.5 * dt, u))
    k3 = np.array(func(state + 0.5 * dt * k2, t + 0.5 * dt, u))
    k4 = np.array(func(state + dt * k3, t + dt, u))

    next_state = state + (dt/6) * (k1 + 2*k2 + 2*k3 + k4)

    return next_state

# 设置仿真参数
dt = 0.1
t = np.arange(0, 20, dt)
t_now = np.zeros_like(t)
x0 = np.zeros_like(t)
y0 = np.zeros_like(t)
z0 = np.zeros_like(t)
q0 = np.zeros_like(t)
q1 = np.zeros_like(t)
q2 = np.zeros_like(t)
q3 = np.zeros_like(t)
theta1 = np.zeros_like(t)
theta2 = np.zeros_like(t)
theta3 = np.zeros_like(t)
theta4 = np.zeros_like(t)
theta5 = np.zeros_like(t)
theta6 = np.zeros_like(t)
vx0 = np.zeros_like(t)
vy0 = np.zeros_like(t)
vz0 = np.zeros_like(t)
wx0 = np.zeros_like(t)
wy0 = np.zeros_like(t)
wz0 = np.zeros_like(t)
dtheta1 = np.zeros_like(t)
dtheta2 = np.zeros_like(t)
dtheta3 = np.zeros_like(t)
dtheta4 = np.zeros_like(t)
dtheta5 = np.zeros_like(t)
dtheta6 = np.zeros_like(t)
control1 = np.zeros_like(t)
control2 = np.zeros_like(t)
control3 = np.zeros_like(t)
control4 = np.zeros_like(t)
control5 = np.zeros_like(t)
control6 = np.zeros_like(t)
control7 = np.zeros_like(t)
control8 = np.zeros_like(t)
control9 = np.zeros_like(t)
control10 = np.zeros_like(t)
control11 = np.zeros_like(t)
control12 = np.zeros_like(t)
t_now[0] = 0.0
x0[0] = 0.0
y0[0] = 0.0
z0[0] = 0.0
q0[0] = 1.0
q1[0] = 0.0
q2[0] = 0.0
q3[0] = 0.0
theta1[0] = 0.0
theta2[0] = 0.0
theta3[0] = 0.0
theta4[0] = 0.0
theta5[0] = 0.0
theta6[0] = 0.0
vx0[0] = 0.0
vy0[0] = 0.0
vz0[0] = 0.0
wx0[0] = 0.0
wy0[0] = 0.0
wz0[0] = 0.0
dtheta1[0] = 0.0
dtheta2[0] = 0.0
dtheta3[0] = 0.0
dtheta4[0] = 0.0
dtheta5[0] = 0.0
dtheta6[0] = 0.0
control1[0] = 0.0
control2[0] = 0.0
control3[0] = 0.0
control4[0] = 0.0
control5[0] = 0.0
control6[0] = 0.0
control7[0] = 0.0
control8[0] = 0.0
control9[0] = 0.0
control10[0] = 0.0
control11[0] = 0.0
control12[0] = 0.0


for i in range(1, len(t)):
    state = np.array([x0[i-1], y0[i-1], z0[i-1], q0[i-1], q1[i-1], q2[i-1], q3[i-1], theta1[i-1], theta2[i-1], theta3[i-1], theta4[i-1], theta5[i-1], theta6[i-1], vx0[i-1], vy0[i-1], vz0[i-1], wx0[i-1], wy0[i-1], wz0[i-1], dtheta1[i-1], dtheta2[i-1], dtheta3[i-1], dtheta4[i-1], dtheta5[i-1], dtheta6[i-1]])
    u = controller(state)
    next_state = runge_kutta_step(system_dynamics, t_now[i-1], state, u, dt)
    t_now[i] = t_now[i-1] + dt
    x0[i], y0[i], z0[i], q0[i], q1[i], q2[i], q3[i], theta1[i], theta2[i], theta3[i], theta4[i], theta5[i], theta6[i], vx0[i], vy0[i], vz0[i], wx0[i], wy0[i], wz0[i], dtheta1[i], dtheta2[i], dtheta3[i], dtheta4[i], dtheta5[i], dtheta6[i] = next_state
    control1[i], control2[i], control3[i], control4[i], control5[i], control6[i], control7[i], control8[i], control9[i], control10[i], control11[i], control12[i]= u

    print(t_now[i-1])


# 绘制结果图形

plt.plot(t, 57.3*wx0, label='wx0')
plt.plot(t, 57.3*wy0, label='wy0')
plt.plot(t, 57.3*wz0, label='wz0')
plt.xlabel('Time')
plt.ylabel('Value')
plt.title('System Response')
plt.legend()
plt.grid(True)
plt.show()

plt.plot(t, vx0, label='vx0')
plt.plot(t, vy0, label='vy0')
plt.plot(t, vz0, label='vz0')
plt.xlabel('Time')
plt.ylabel('Value')
plt.title('System Response')
plt.legend()
plt.grid(True)
plt.show()

plt.plot(t, q0, label='q0')
plt.plot(t, q1, label='q1')
plt.plot(t, q2, label='q2')
plt.plot(t, q3, label='q3')
plt.xlabel('Time')
plt.ylabel('Value')
plt.title('System Response')
plt.legend()
plt.grid(True)
plt.show()

plt.plot(t, x0, label='x0')
plt.plot(t, y0, label='y0')
plt.plot(t, z0, label='z0')
plt.xlabel('Time')
plt.ylabel('Value')
plt.title('System Response')
plt.legend()
plt.grid(True)
plt.show()

plt.plot(t, 57.3*dtheta1, label='dtheta1')
plt.plot(t, 57.3*dtheta2, label='dtheta2')
plt.plot(t, 57.3*dtheta3, label='dtheta3')
plt.plot(t, 57.3*dtheta4, label='dtheta4')
plt.plot(t, 57.3*dtheta5, label='dtheta5')
plt.plot(t, 57.3*dtheta6, label='dtheta6')
plt.xlabel('Time')
plt.ylabel('Value')
plt.title('System Response')
plt.legend()
plt.grid(True)
plt.show()
# 绘制结果图形
plt.plot(t, control1, label='T1')
plt.plot(t, control2, label='T2')
plt.plot(t, control3, label='T3')
plt.plot(t, control4, label='T4')
plt.plot(t, control5, label='T5')
plt.plot(t, control6, label='T6')
plt.plot(t, control7, label='T7')
plt.plot(t, control8, label='T8')
plt.plot(t, control9, label='T9')
plt.plot(t, control10, label='T10')
plt.plot(t, control11, label='T11')
plt.plot(t, control12, label='T12')
plt.xlabel('Time')
plt.ylabel('T')
plt.title('System Response')
plt.legend()
plt.grid(True)
plt.show()