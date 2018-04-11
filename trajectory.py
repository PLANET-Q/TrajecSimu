#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Sat Jul  8 23:06:17 2017

@author: shugo
"""

import numpy as np
from scipy.integrate import ode, odeint
from scipy import fftpack, interpolate, integrate, optimize
# import pandas as pd
import quaternion 
from Rocket_simu import Rocket
import matplotlib.pyplot as plt

class trajec_main(Rocket):
    """
    ====================================================
    This class is a sub-class for main ODE trajecotry computation
    ====================================================
    """
    
    """
    ----------------------------------------------------
        Method for initial setup    
    ----------------------------------------------------
    """
    def __init__(self, params_df):
        # =============================================
        # this method is called when instance is created. Parameter setup is also done by this method.
        #
        # INPUT: params_df = parameters contains all parameters needed
        # =============================================
        # setup parameters in the instance by calling superclasses method
        self.get_default()   
        self.overwrite_parameters(params_df)   
        
        # initialize a dict for result
        self.res_trajec_main = {'v_para_deploy': -1, 'v_launch_clear': -1}

        return None     

        
    """
    ----------------------------------------------------
        Method for main ODE integration  
    ----------------------------------------------------
    """
    
    def ODE_main(self):
        # =============================================
        # this method runs ODE integration given all the parameters
        #     
        # OUTPUT: solution: shape (len(t), len(y0))
        #                   containing the value of y for each desired time in t, with the initial value y0 in the first row.
        # =============================================
        
        # ---------------------------------------------
        #      Initialization
        # ---------------------------------------------
        # set initial value
        t0 = 0.              # start time
        u0 = self.SetupICs() # initiate state vector

        # number of iteration
        self.N_iter = 0
            
        # initialize a list for backup
        # backup = [time,flag,u]
        self.backup = [np.r_[t0,0,u0]]

        # landing time initialization: 
        self.landing_time = self.t_max
            
        # set flag = 1 (on launch rail)
        self.flag = 1
        
        # setup Cd0 interpolation curve
        self.setup_Cd0()
    
        """
        print('----------------------------')
        print('  We are GO for launch')
        print('----------------------------')
        print(' ')
        """
          
        """
        self.record_tmp = np.array([0]) 
        self.time_all = np.array([0]) 
        """
        
        if self.integ == 'lsoda_odeint':
            # ---------------------------------------------
            #      use scipy.odeint
            # ---------------------------------------------
            # create time array
            
            self.t = np.arange(t0,self.t_max,self.dt)
            #try:
                # ---  run trajectory computation   ---
            self.solution = odeint(self.f_main, u0, self.t)
            #except:
            #    print('ode integration failed!')
            #    pass
            
        else:
            
            print('currently not supported. set "inter,lsoda_odeint" in csv file.')
            
            """
            # ---------------------------------------------
            #      use scipy.ode
            # ---------------------------------------------            
            # create ODE integration instance of scipy.integrate.ode
            r = ode(self.f_main).set_integrator(self.integ)
            # set initial value
            r.set_initial_value(u0, t0)
            
            # ---  run trajectory computation   ---
            # loop until landing (or maximum time) 
            while r.successful() and r.t < self.t_max and self.flag<5:
                # call ode integration
                r.integrate(r.t+self.dt)
                #END IF
            #END WHILE
            
            """
        #END IF
        
        """        
        print('----------------------------')
        print('   Completed trajectory computation ')
        print('----------------------------')
        print(' ')
        print(' ')
        """ 
        """
        plt.figure
        plt.plot(self.time_all, self.record_tmp)
        # plt.plot(self.t_MECO, self.m_dry, '*')
        plt.grid()
        plt.show()
        """
        
        return None   
        

    """
    ----------------------------------------------------
        Other methods for simulation      
    ----------------------------------------------------
    """
    
    def SetupICs(self):
        # =======================================
        # this method sets up initial conditions
        # and returns initial state vector u0: 13*1 array
        #  u0[0:3]  = (1*3) vector: translation xyz      referring to local fixed coord.
        #  u0[3:6]  = (1*3) vector: velocity uvw         referring to body coord. 
        #  u0[6:10] = quaternion:   attitude             convert from local to body
        #  u0[10:]  = (1*3) vector: abgular velocity pqr referring to body coord. 
        # =======================================
        
        # initial location
        x0 = np.zeros(3)
         
        # initial velocity
        v0 = np.zeros(3)
        
        # initial angular velocity
        omega0 = np.zeros(3)
        
        # compute initial quaternion 
        # convert body coord. to local fixed coord.
        angle1 = (90.-self.azimuth)/2. * (np.pi/180.)
        qz = np.quaternion(np.cos(angle1),0.,0.,np.sin(angle1))  # around z
        angle2 = -self.elev_angle/2. * (np.pi/180.)
        qy = np.quaternion(np.cos(angle2),0.,np.sin(angle2),0.)  # around y
        # convert local -> rotate around z -> rotate around y -> body
        q0 = qz * qy
        
        # convert quaternion to float array
        q01 = quaternion.as_float_array(q0)
        
        # initial state "vector"
        u0 = np.r_[x0,v0,q01,omega0]
    
        return u0
        

    def f_main(self,u,t):
        # =======================================
        # this method is the Right hand side function of ODE (du_dt = f_main(u))
        #
        # INPUT: t    = time (scalar) 
        #        u    = state vector u (13*1)
        #
        # NOTE:  self.flag = 0 before ignition
        #                    1   : on launch rail
        #                    2   : thrusted flight
        #                    3   : inertial flight
        #                    3.5 : drogue chute deployed, before main para deployment
        #                    4   : main parachute deployed
        #                    5   : landed
        # =======================================
        
        # if the rocket has already landed, return 0
        if self.flag == 5:
            return u*0.
            
        # swap input when we use scipy.integrate.ode
        if self.integ != 'lsoda_odeint':
            tmp = u
            u = t
            t = tmp
        #END IF
    
        # count the number of function call
        self.N_iter += 1
        
        # backup
        if np.mod(self.N_iter, self.N_record) == 0:
            self.add_backup(t,u)
        #END IF
        
        # --------------------------
        #   extract vectors
        # --------------------------
        # x =     translation         :expressed in local-fixed coordinate
        # v =     velocity            :expressed in body coordinate
        # q =     atitude quaternion  :conversion from local-fixed to body
        # omega = angular velocity    :expressed in body coordinate
        x,v,q,omega = self.state2vecs_quat(u)
         
            
        # ----------------------------
        #    Direction Cosine Matrix for input q
        # ----------------------------
        # Tbl = transform from local(fixed) coordinate to body coord.
        #     note: as_rotation_matrix is for vector rotation
        #         -> for coordinate rotation, input conj(q)
        Tbl = quaternion.as_rotation_matrix(np.conj(q))
        
        # ----------------------------
        #   flight mode classification
        # ----------------------------
        if self.flag==1 and x[2] > self.rail_height:
            # detect launch clear

            # get air speed and AoA
            air_speed, _, AoA, _ = self.air_state(x,v,Tbl)
            print('----------------------------')
            print('  Launcher-clear at t = ',np.round(t, 2),'[s]')
            print('  ground speed: ',np.round(np.linalg.norm(v), 2),'[m/s]')
            print('  true air speed: ', np.round(air_speed, 2), '[m/s], AoA: ', np.round(AoA*180./np.pi, 3), '[deg]')
            print('----------------------------')
            print(' ')
            
            # record history
            self.add_backup(t,u)
            # record launch clear air speed
            self.res_trajec_main.update( { 'v_launch_clear' : air_speed } )
            # switch into thrusted flight
            self.flag = 2   
            
        elif (self.flag==1 or self.flag==2) and t >= self.t_MECO:
            # detect MECO

            print('----------------------------')
            print('  MECO at t = ',np.round(t, 2), '[s]')
            print('  current altitude: ',np.round(x[2], 2), '[m]')
            print('  ground speed:    ',np.round(np.linalg.norm(v), 2), '[m/s]')
            print('----------------------------')
            print(' ')
            
            # record history
            self.add_backup(t,u)
            # switch into coasting flight
            self.flag = 3   
            
        elif self.flag==3 and t >= self.t_deploy:
            # detect 1st parachute deployment

            # get air speed and AoA
            air_speed, _, AoA, _ = self.air_state(x,v,Tbl)
            print('----------------------------')
            print('  1st parachute deployed at t = ', np.round(t, 2), '[s]')
            print('  current altitude: ', np.round(x[2], 2), '[m]')
            print('  ground speed:    ', np.round(np.linalg.norm(v), 2), '[m/s]')
            print('  true air speed: ', np.round(air_speed, 2), '[m/s], AoA: ', np.round(AoA*180./np.pi, 3), '[deg]')
            print('----------------------------')
            print(' ')
            
            # record history
            self.add_backup(t,u)
            # record parachute deploy air speed
            self.res_trajec_main.update( {'v_para_deploy' : air_speed } )
            
            # switch into parachute fall
            if self.flag_2ndpara:
                # when two stage separation is True
                self.flag = 3.5  #flag: drogue chute deployed
            else:
                self.flag = 4   # flag: main chute deployed
            
            # stop rotation
            omega = np.zeros(3)  
            
        # elif self.flag==3.5 and (t >= self.t_deploy_2 or 高度が規定値以下) :
        elif self.flag==3.5 and t >= self.t_deploy_2:
            # detect 2nd parachute deployment

            # get air speed and AoA
            air_speed, _, AoA, _ = self.air_state(x,v,Tbl)
            print('----------------------------')
            print('  2nd parachute deployed at t = ', np.round(t, 2), '[s]')
            print('  current altitude: ', np.round(x[2], 2), '[m]')
            print('  ground speed:    ', np.round(np.linalg.norm(v), 2), '[m/s]')
            print('  true air speed: ', np.round(air_speed, 2), '[m/s]' )
            print('----------------------------')
            print(' ')
            
            # overwrite parachute properties with 2nd para
            self.Cd_para = self.Cd_para_2
            self.S_para = self.S_para_2
            
            # record history
            self.add_backup(t,u)
            
            # switch into main parachute fall
            self.flag = 4
        
            
        elif self.flag > 1 and self.flag < 5 and x[2] < 0. :
            # detect landing
            print('----------------------------')
            print('  Landing at t = ',np.round(t, 2), '[s]')
            print('  landing ground speed: ', np.round(np.linalg.norm(v), 2), '[m/s]')
            print('          location x = ', np.round(x[0], 2), '[m]')
            print('                   y = ', np.round(x[1], 2), '[m]')
            print('----------------------------')
            print(' ')
            # record landing time
            self.landing_time = t
            self.res_trajec_main.update( {'landing_time' : t } )
            
            # record history
            self.add_backup(t,u)
            # flag -> landed
            self.flag = 5 
            # quit integration
            if self.integ == 'lsoda_odeint':
                return u*0.                
        #END IF 
        
        
        # ----------------------------
        #    1. Translation
        # ----------------------------
        # translation time rate = velocity
        # convert v from body coord. to fixed coord.
        dx_dt = np.dot(Tbl.T,v)
        
        # call mass properties 
        mass,MOI,d_dt_MOI,CG = self.mass_MOI(t)  # note that jet-damping effect is included in d(MOT)/dt
        
        # ----------------------------
        #    2. Velocity
        # ----------------------------
        # velocity time rate = mass * force
        # force = rotation effect + aero + thrust + gravity 
        #      where  aero   = aero(u):   function of state u
        #             thrust = thrust(t): function of time t
        #             mass   = mass(t):   function of time t   
    
        # get aerodynamic force/moment
        aeroF, aeroM = self.aero(x,v,omega,Tbl,CG)  
        # self.time_all = np.append(self.time_all, t)
        # self.record_tmp = np.append(self.record_tmp, aeroF[0])
    
        # set external force depending on the state
        if self.flag == 1:
            # on launch rail -> du/dx only. Consider rail-rocket friction force 
            # total acceleration 
            dv_dt = -np.cross(omega,v) + np.dot(Tbl, self.grav) + (aeroF + self.thrust(t) + self.friction() ) / mass
            # cancell out y,z
            dv_dt = np.array([dv_dt[0],0.,0.])
            
            if dv_dt[0] < 0:
                # when du/dx is negative (i.e. weight is greater that thrust)
                # -> rocket is hold up. return zero acceleration
                dv_dt = np.zeros(3)
            #END
            
        elif self.flag == 2:
            # thrust ON
            # total acceleration 
            dv_dt = -np.cross(omega,v) + np.dot(Tbl, self.grav) + self.Coriolis(v, Tbl) + (aeroF + self.thrust(t)) / mass
            
        elif self.flag == 3 or self.flag == 5:
            # coasting phase
            # total acceleration 
            dv_dt = -np.cross(omega,v) + np.dot(Tbl, self.grav) + self.Coriolis(v, Tbl) + aeroF / mass

        elif self.flag == 3.5 or self.flag == 4:
            # parachute deployed
            dv_dt = np.dot(Tbl, self.grav) + self.Coriolis(v, Tbl) + self.parachute_F(x,v,Tbl) / mass
        #END IF
        
        
        # ----------------------------
        #    3. Atitude 
        # ----------------------------
        # represented in quaternion which rotates from fixed to body
        # quaternion differential equation
        
        # stop rotation when parachute ihas deployed
        if self.flag == 3.5 or self.flag == 4:
            omega *= 0.
        #END IF
        
        # convert omega to quaternion
        q_omega = np.r_[[0.], omega]
        q_omega2 = quaternion.as_quat_array(q_omega)
        
        # dq/dt to be returned
        dq_dt2 = 0.5 * q * q_omega2
        
        # convert back to float array
        dq_dt = quaternion.as_float_array(dq_dt2)
        
        # ----------------------------
        #    4. Angular velocity
        # ----------------------------
        # angular velocity time rate: comes from Euler equation
        #      where  MOI       = MOI(t):       function of time
        #             d_dt(MOI) = d_dt(MOI)(t): function of time
        #             aeroM     = aeroM(u):     function of state 
        #             ThrustM   = ThrustM(t)=0: function of time
        
        if self.flag == 1 or self.flag == 3.5 or self.flag == 4:
            # on launch rail /parachute deployed -> no angular velocity change
            domega_dt = np.zeros(3)

        else:
            # MOI1 = MOI(t)           # moment of inertia     
            # d_dt_MOI1 = d_dt_MOI(t) # time derivative of MOI
            # aeroM1 = aeroM(u)       # aerodynamic moment
            
            # Euler eqn of rotation
            #tmp1 = 1./MOI[0] * ( (MOI[1]-MOI[2])*omega[1]*omega[2] - d_dt_MOI[0]*omega[0] + aeroM[0]) 
            #tmp2 = 1./MOI[1] * ( (MOI[2]-MOI[0])*omega[0]*omega[2] - d_dt_MOI[1]*omega[1] + aeroM[1]) 
            #tmp3 = 1./MOI[2] * ( (MOI[0]-MOI[1])*omega[0]*omega[1] - d_dt_MOI[2]*omega[2] + aeroM[2]) 
            #domega_dt = np.array([tmp1,tmp2,tmp3])
            domega_dt = (-np.cross(omega,MOI*omega) - d_dt_MOI*omega + aeroM) / MOI
        # END IF

        # ----------------------------
        #    Set variables back in the state vector form
        # ----------------------------
        du_dt = np.r_[dx_dt,dv_dt,dq_dt,domega_dt]

        # ----------------------------
        #      apogee detection 
        # ----------------------------
        if dx_dt[2] < 0.:
            self.apogee_count += 1
            
            if self.apogee_count >= 10:
                # detect apogee and set parachute deployment time
                # parachute is deployed at t = t_deploy, 
                # which is either "t_deploy" [s] after ignition or
                # "t_para_delay" [s] after apogee detection
                
                self.t_deploy = min( self.t_deploy, t + self.t_para_delay )
                self.apogee_count = -1.e10
            #END IF
        #END IF
        
        return du_dt
        
        
    def mass_MOI(self,t):
        # =======================================
        # this method returns mass properties of rocket
        # 
        # INPUT:  t = time
        # OUTPUT: mass = total mass of rocket
        #         MOI = total moment of inertia wrt CG
        #         d_dt_MOI = MOI time rate
        #         CG = center of gravity location from the nose tip
        # =======================================

        if t >= self.t_MECO:
            # ---------------------------
            # mass for coasting phase (m, I = const.)
            # ---------------------------
            # mass
            mass = self.m_dry
            # moment of inertia
            MOI = self.MOI_dry
            # total CG location
            CG = self.CG_dry
            
            return mass, MOI, np.zeros(3), CG
            
        else:
            # ---------------------------
            # mass for powered phase (m = m(t), I = I(t))
            # ---------------------------
            # propellant comsumption rate = (impulse consumed so far) / (total impulse)
            time_so_far = np.linspace(0., t)
            Impulse_so_far = integrate.trapz(self.thrust_function(time_so_far), time_so_far)
            r = ( 1 - Impulse_so_far/self.Impulse_total )  # impulse ratio
            
            # total mass
            mass = self.m_dry + r * self.m_prop
            # total CG location
            CG = (self.CG_dry*self.m_dry + self.CG_prop*r*self.m_prop) / mass
            # total MOI using parallel axis theorem
            tmp = np.array([0.,1.,1.])
            MOI = self.MOI_dry + tmp*self.m_dry*(CG-self.CG_dry)**2. + r*self.MOI_prop + tmp*(CG-self.CG_prop)*(r*self.m_prop)**2.
    
            # ---------------------------------
            # finite differencing for d(MOI)/dt, dm/dt
            # ---------------------------------
            h = 1.E-3
            Impulse_so_far = integrate.trapz(self.thrust_function(np.linspace(0., t+h)), np.linspace(0., t+h))
            r2 = (1- Impulse_so_far/self.Impulse_total)  # impulse ratio
            # total mass
            # mass2 = self.m_dry + r2 * self.m_prop
            # total CG location
            CG2 = (self.CG_dry*self.m_dry + self.CG_prop*r2*self.m_prop) / (self.m_dry + r2*self.m_prop) 
            # total MOI using parallel axis theorem
            MOI2 = self.MOI_dry + tmp*self.m_dry*(CG2-self.CG_dry)**2. + r2*self.MOI_prop + tmp*(CG2-self.CG_prop)*(r2*self.m_prop)**2.

            Impulse_so_far = integrate.trapz(self.thrust_function(np.linspace(0., t-h)), np.linspace(0., t-h))
            r3 = (1- Impulse_so_far/self.Impulse_total)  # impulse ratio
            # mass3 = self.m_dry + r3 * self.m_prop
            CG3 = (self.CG_dry*self.m_dry + self.CG_prop*r3*self.m_prop) / (self.m_dry + r3*self.m_prop)
            MOI3 = self.MOI_dry + tmp*self.m_dry*(CG3-self.CG_dry)**2. + r3*self.MOI_prop + tmp*(CG3-self.CG_prop)*(r3*self.m_prop)**2.
    
            # dm/dt and d(MOI)/dt
            # d_dt_m = (mass2 - mass3) / (2*h)
            d_dt_MOI = (MOI2 - MOI3) / (2*h)
            
            return mass, MOI, d_dt_MOI, CG
        #END IF
                
    def Coriolis(self, v, Tbl):
        # =======================================
        # computes Coriolis force per unit mass
        # 
        # INPUT:  v = velocity in body coord.
        #         Tbl = tranformation matrix from local coord. to body coord.
        # OUTPUT: Fcor = Coriolis force vector in body coord.
        # =======================================
        
        # Coriolis  force in body coord. note that self.omega_earth, omega of earth-spin, is given in local coord.
        Fcor = 2. * np.cross( (Tbl@self.omega_earth), v )
        
        return Fcor
    
    def thrust(self,t):
        # =======================================
        # returns thrust force 
        # 
        # INPUT:  t = time
        # OUTPUT: T = thrust vector in body coord.
        # =======================================
            
        # get thrust power at t=t
        tmp = self.thrust_function(t)
        if tmp < 0:
            # if thrust is nagative, which might happen because of curve fitting, overwrite with 0
            tmp = 0.
            
        # thrust vector (assume 0 misalignment)
        T = np.array([tmp,0.,0.])
        
        return T
        
        
    def friction(self):
        # =======================================
        # returns friction force (rail-rocket contact) 
        # 
        # OUTPUT: friction = friction force vector in body coord.
        # =======================================
        friction = np.zeros(3)
        
        return friction
            
        
            
    def aero(self,x,v,omega,Tbl,CG):
        # =======================================
        # returns aerodynamic force and moment 
        #
        # INPUT:  x     = translation         :in fixed coordinate
        #         v     = velocity            :in body coordinate
        #         omega = angular velocity    :in body coordinate
        #         Tbl   = transform matrix from local(fixed) coordinate to body coord.
        #         CG    = current CG location
        # OUTOUT: force_all  = aerodynamic force vector in body coord.
        #         moment_all = aerodynamic moment vector wrt. CG in body coord.
        # =======================================
        
        # --------------------------------------------------
        #   Compute air velocity, angle-of-attack, roll-angle
        # --------------------------------------------------
        air_speed, v_air, alpha, phi = self.air_state(x,v,Tbl)
            
        # ----------------------------------
        #   aerodynamic force on body ( might exclude fins)
        # ----------------------------------
        # air property at the altitude
        _,_,rho,a = self.standard_air(x[2])        
        
        # Mach number
        Mach = air_speed/a
    
        # drag/lift coefficient
        Cd,Cl = self.rocket_coeff_nofins(Mach,alpha)
        
        # convert coefficient to body coord.
        cosa = np.cos(alpha)
        sina = np.sin(alpha)
        
        #C1b = -Cl*sina + Cd*cosa
        #C2b = -(Cl*cosa + Cd*sina)*np.sin(phi)
        #C3b = (Cl*cosa + Cd*sina)*np.cos(phi)
        C = np.array([ (-Cl*sina + Cd*cosa), \
                      (Cl*cosa + Cd*sina)*np.sin(phi), \
                      (Cl*cosa + Cd*sina)*np.cos(phi)])    # need check here
        
        # force acting on CP of the body
        force_all = 0.5 * rho * air_speed**2. * self.X_area * (-C)
        
        # aerodynamic moment wrt CG
        moment_all = np.cross( np.array([CG-self.CP_body,0.,0.]) , force_all )
        
        # add aerodynamic damping effect  (note: Cm_omega < 0)
        moment_all += 0.25 * rho * air_speed * self.Cm_omega_bar * omega
        
        if self.aero_fin_mode == 'indiv':
            # ----------------------------------
            #   compute force on fins individually
            # ----------------------------------
            # force and moment (wrt fin Leading edge)
            force_fin, moment_fin_tmp = self.fin_aero(v_air,omega[0],rho)
            # moment generated by fin wrt CG
            moment_fin = moment_fin_tmp + np.cross( np.array([CG-self.LE_fins,0.,0.]) , force_fin )
            
            # total force
            force_all += force_fin
            # total moment
            moment_all += moment_fin
        # END IF 
        
        """
        F_fill_1d = 0.5 * rho * air_speed**2. * 0.03*0.1 * 0.7
        F_fill = np.array([-F_fill_1d, 0., 0.])
        F_fill_mom = np.cross( np.array([0.,0.13, 0.]), F_fill)
        
        force_all += F_fill
        moment_all += F_fill_mom
        """
        
        return force_all, moment_all
        
    def fin_aero(self,v_air,omega_x,rho):
        # ==============================================
        # return aerodynamic force and moment generated by fins
        # reference point of moment = [Leading Edge * root] point
        # ==============================================
        
        if all(v_air==0):
            v_air = np.array([0.0001,0.,0.])
            
        # each component of air velocity 
        u = v_air[0]  
        v = v_air[1]
        w = v_air[2]
        
        # --------------------
        #    fin1: z=0, y>0
        # --------------------
        # net air speed on xz plane. (vertical to the fin) Take roll-spin into account.
        U1 = np.sqrt(u**2. + ( w - self.r_arm*omega_x )**2.)  # array as func. of r_arm
                    
        # net angle-of-attack (air flow, roll-spin, fin attachment angle are considered)
        alpha1 = np.arcsin(( w - self.r_arm*omega_x )/U1) + self.alpha_fins
    
        # lift force distribution
        lift1 = np.pi * rho * U1**2. * self.fin_len * np.sin(alpha1) * self.dy_fin  
        # convert into body coordinate
        force1x = - lift1 * np.sin(alpha1)
        force1z = lift1 * np.cos(alpha1)
        
        # roll moment distribution
        rollmoment1 = force1z * self.r_arm
        # pitch moment distribution
        pitchmoment1 = force1z * self.xFP
        
        # total force of fin1
        f1 = np.array([np.sum(force1x),0.,np.sum(force1z)])
        
        # total moment of fin1
        m1 = np.array([np.sum(rollmoment1),np.sum(pitchmoment1),0.])
        
        
        # --------------------
        #    fin2: y=0, z>0
        # --------------------
        # net air speed on xz plane. (vertical to the fin) Take roll-spin into account.
        U2 = np.sqrt(u**2. + ( v + self.r_arm*omega_x )**2.)  # array as func. of r_arm
        
        # net angle-of-attack (air flow, roll-spin, fin attachment angle are considered)
        alpha2 = np.arcsin(( v + self.r_arm*omega_x )/U2) - self.alpha_fins
    
        # lift force distribution
        lift2 = np.pi * rho * U2**2. * self.fin_len * np.sin(alpha2) * self.dy_fin  
        # convert into body coordinate
        force2x = - lift2 * np.sin(alpha2)
        force2y = lift2 * np.cos(alpha2)
        
        # roll moment distribution
        rollmoment2 = - force2y * self.r_arm  # positive y-force generates negative roll moment
        # pitch moment distribution
        yawmoment2 = - force2y * self.xFP  # positive y-force generates negative yaw moment
        
        # total force of fin1
        f2 = np.array([np.sum(force2x),np.sum(force2y),0.])
        
        # total moment of fin1
        m2 = np.array([np.sum(rollmoment2),0.,np.sum(yawmoment2)])
        
    
        # --------------------
        #    fin3: z=0, y<0
        # --------------------
        # net air speed on xz plane. (vertical to the fin) Take roll-spin into account.
        U3 = np.sqrt(u**2. + ( w + self.r_arm*omega_x )**2.)  # array as func. of r_arm
        
        # net angle-of-attack (air flow, roll-spin, fin attachment angle are considered)
        alpha3 = np.arcsin(( w + self.r_arm*omega_x )/U3) - self.alpha_fins
    
        # lift force distribution
        lift3 = np.pi * rho * U3**2. * self.fin_len * np.sin(alpha3) * self.dy_fin  
        # convert into body coordinate
        force3x = - lift3 * np.sin(alpha3)
        force3z = lift3 * np.cos(alpha3)
        
        # roll moment distribution
        rollmoment3 = -force3z * self.r_arm  # positive z-force generates negative moment
        # pitch moment distribution
        pitchmoment3 = force3z * self.xFP
        
        # total force of fin3
        f3 = np.array([np.sum(force3x),0.,np.sum(force3z)])
        
        # total moment of fin3
        m3 = np.array([np.sum(rollmoment3),np.sum(pitchmoment3),0.])
        
        # --------------------
        #    fin4: y=0, z<0
        # --------------------
        # net air speed on xz plane. (vertical to the fin) Take roll-spin into account.
        U4 = np.sqrt(u**2. + ( v - self.r_arm*omega_x )**2.)  # array as func. of r_arm
        
        # net angle-of-attack (air flow, roll-spin, fin attachment angle are considered)
        alpha4 = np.arcsin(( v - self.r_arm*omega_x )/U2) + self.alpha_fins
    
        # lift force distribution
        lift4 = np.pi * rho * U4**2. * self.fin_len * np.sin(alpha4) * self.dy_fin  
        # convert into body coordinate
        force4x = - lift4 * np.sin(alpha2)
        force4y = lift4 * np.cos(alpha2)
        
        # roll moment distribution
        rollmoment4 = force4y * self.r_arm  # positive y-force generates positive roll moment
        # pitch moment distribution
        yawmoment4 = - force4y * self.xFP  # positive y-force generates negative yaw moment
        
        # total force of fin1
        f4 = np.array([np.sum(force4x),np.sum(force4y),0.])
        
        # total moment of fin1
        m4 = np.array([np.sum(rollmoment4),0.,np.sum(yawmoment4)])
        
        # ---------------------
        #   total force/moment by fins
        # ---------------------
        force_all = f1+f2+f3+f4
        moment_all = m1+m2+m3+m4
        
        return force_all, moment_all

    def air_state(self,x,v,Tbl):
        # =======================================
        # returns air speed, air velocity, angle of attack, and roll angle
        #
        # INPUT:  x            = translation         :in fixed coordinate
        #         v            = velocity            :in body coordinate
        #         Tbl          = transform matrix from local(fixed) coordinate to body coord.
        # OUTOUT: air_speed    = air speed [m/s]
        #         v_air        = air velocity vector
        # =======================================

        # air velocity = -rocket_velocity + wind_velocity
        # NOTE: wind(x) is wind velocity in local coord. need conversion to body coord.
        v_air = -v + np.dot( Tbl, self.wind(x[2]) )   
        air_speed = np.linalg.norm(v_air) # air speed (positive scalar)

        # total angle-of-attack
        if air_speed == 0:
            alpha = 0.
        else:
            alpha = np.arccos( -v_air[0]/air_speed )
            #if v_air[0] > 0:
            #    alpha = -alpha
            #END IF
        #END IF
        
        # if np.isnan(alpha):
        #    alpha = 0.

        # roll-angle
        if v_air[2]==0:
            # if w = 0, atan(v/w) is not defined 
            if -v_air[1] > 0:
                phi = np.pi /2.
            else:
                phi = -np.pi /2.    
        else:
            phi = np.arctan( -v_air[1]/ -v_air[2] )
            if -v_air[2]<0:
                # if w<0, that is in 3rd or 4th quadrant in yz-plane, add pi to phi
                #  so that sign of sin(phi), cos(phi) will be correctly calculated.
                phi += np.pi
            #END IF
        #END IF
                
        #if np.isnan(phi):
        #    phi = np.arctan( v_air[1]+0.0001/v_air[2]+0.0001 )

        return air_speed, v_air, alpha, phi
    
    
    def parachute_F(self,x,v,Tbl):
        # =======================================
        # returns aerodynamic force and moment 
        #
        # INPUT:  x     = translation         :in fixed coordinate
        #         v     = velocity            :in body coordinate
        #         Tbl   = transform matrix from local(fixed) coordinate to body coord.
        # OUTOUT: parachute_drag  = parachute drag force in body coord.
        # =======================================
        
        # air velocity = -rocket_velocity + wind_velocity
        #    NOTE: wind(x) is wind velocity in local coord. need conversion to body coord.
        v_air = -v + np.dot( Tbl,self.wind(x[2]) )   

        # air property at the altitude
        _,_,rho,_ = self.standard_air(x[2])

        # parachute drag force 
        parachute_drag = 0.5 * rho * np.linalg.norm(v_air) * v_air * self.S_para * self.Cd_para

        # print('v_air', v_air, 'speed', np.linalg.norm(v_air))
        return parachute_drag

        
    def standard_air(self,h):
        # ==============================================
        # returns air property given an altitude 
        # INPUT: h = altitude [m]
        # ==============================================
        
        # gas constant
        R = 287.15  # [J/kg.K]
        # gravitational accel.
        g = 9.81  # [m/s^2]
        
        
        """
        # temperature goes down 0.0065K/m until it reaches -56.5C (216.5K)
        #                                       it is approximately 11km
        T = self.T0 - 0.0065*h # [K]
        
        
        # temperature is const at 216.5 K for alt. < 20km
        if type(T) == np.ndarray:
            T[T<216.5] = 216.5
            
        elif T < 216.5:
            T = 216.5
        
        # pressure
        p = self.p0 * (T/self.T0)**5.256  #[Pa]
        """
        
        
        if h <= 11.*10**3:
            # *** Troposphere ***
            # temperature lapse rate
            gamma = -0.0065
            # temperature 
            T = self.T0 + gamma * h # [K]
            #pressure
            p = self.p0 * (T/self.T0)**(-g/(gamma*R)) #[Pa]
            
        elif h <= 20.*10**3:
            # *** Tropopause ***
            # temperature is const at 11km-20km
            # p11 = pressure at 11km alt.
            T,p11,_,_ = self.standard_air(11000.)
            # pressure 
            p = p11 * np.exp( (-g/(R*T)) * (h-11000.) )
            
        elif h <= 32.*10**3:
            # *** Stratosphere 1 ***
            # temp, pressure at 20km alt.
            T20,p20,_,_ = self.standard_air(20000.)
            # temperature lapse rate
            gamma = 0.001
            # temperature 
            T = T20 + gamma * (h-20000.) # [K]
            #pressure
            p = p20 * (T/T20)**(-g/(gamma*R)) #[Pa] 
            
        elif h <= 47.*10**3:
            # *** Stratosphere 2 ***
            # temp, pressure at 32km alt.
            T32,p32,_,_ = self.standard_air(32000.)
            # temperature lapse rate
            gamma = 0.0028
            # temperature 
            T = T32 + gamma * (h-32000.) # [K]
            #pressure
            p = p32 * (T/T32)**(-g/(gamma*R)) #[Pa] 
            
        elif h <= 51.*10**3:
            # *** Stratopause ***
            # temp, pressure at 47km alt.
            T,p47,_,_ = self.standard_air(47000.)
            # pressure 
            p = p47 * np.exp( (-g/(R*T)) * (h-47000.) )
        
        elif h <= 71.*10**3:
            # *** Mesosphere 1 ***
            # temp, pressure at 51km alt.
            T51,p51,_,_ = self.standard_air(51000.)
            # temperature lapse rate
            gamma = -0.0028
            # temperature 
            T = T51 + gamma * (h-51000.) # [K]
            #pressure
            p = p51 * (T/T51)**(-g/(gamma*R)) #[Pa] 
            
        elif h <= 85.*10**3:
            # *** Mesosphere 2 ***
            # temp, pressure at 51km alt.
            T71,p71,_,_ = self.standard_air(71000.)
            # temperature lapse rate
            gamma = -0.002
            # temperature 
            T = T71 + gamma * (h-71000.) # [K]
            #pressure
            p = p71 * (T/T71)**(-g/(gamma*R)) #[Pa] 
            
        elif h <= 90.*10**3:
            # *** Mesopause ***
            # temp, pressure at 47km alt.
            T,p85,_,_ = self.standard_air(85000.)
            # pressure 
            p = p85 * np.exp( (-g/(R*T)) * (h-85000.) )
            
        elif h <= 110.*10**3:
            # *** Thermosphere  ***
            # temp, pressure at 51km alt.
            T90,p90,_,_ = self.standard_air(90000.)
            # temperature lapse rate
            gamma = 0.0026675
            # temperature 
            T = T90 + gamma * (h-90000.) # [K]
            #pressure
            p = p90 * (T/T90)**(-g/(gamma*R)) #[Pa] 
            
        else:
            T110,p110,_,_ = self.standard_air(110000.)
            
        #END IF
        
            
        # density
        rho = p/(R*T) #[kg/m^3]
        
        # acoustic speed
        a = np.sqrt(1.4*R*T) # [m/s]
        
        return T,p,rho,a 
        
    
    def wind(self,h):
        # ==============================================
        # returns wind vector given an altitude
        #  follows "wind profile power law"
        #
        # INPUT: h = altitude [m]
        # ==============================================
        
        # wind velocity in local fixed coordinate
        wind_vec = self.wind_unitvec * self.wind_speed * (h/self.wind_alt_std)**self.Cwind  

        return wind_vec
        
        
        
    def rocket_coeff_nofins(self, Mach,alpha):
        # input: Mach number
        #        alpha: angle of attack[rad]
        # -------------------
        # Lift Coefficient
        # -------------------
        
        if self.aero_fin_mode == 'indiv':
            # lift coefficient slope for fin-body individual computations
            #  : slope for body Cl
            k1 = 1.
        # END IF
        Cl = self.Cl_alpha * np.sin(2.*alpha)
        
        
        # -------------------
        # Drag Coefficient
        # -------------------
        # get zero angle Cd
        Cd0 = self.f_cd0(Mach)
        
        if self.flag <= 2:
            # powered phase Cd0 correction
            Cd0 -= 0.1 
        # END IF
        
        Cd0 *= self.Cd0/0.7
        
        # drag coefficient "amplitude" for cosign curve fit
        Cd_bar = 15.
        Cd = Cd0 + Cd_bar*( np.cos(2*alpha + np.pi) +1.)
        
        return Cd, Cl
        

    def setup_Cd0(self):
        # setup Cd0 curve as a function of Mach number
        
        # array from CFD
        Cd0_array = np.array([0.6, 0.706,0.699,0.707,0.756,0.964,1.120,1.094,1.060,1.023,0.985,0.951,0.914,0.877,0.839,0.802,0.7])
        # Cd0_array = np.array([0.7,0.806,0.799,0.810,0.869,1.121,1.312,1.280,1.240,1.195,1.150,1.108,1.065,1.020,0.975,0.930, 0.8])
        Mach_array = np.array([0., 0.6, 0.7, 0.8, 0.9, 1.0, 1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.7, 1.8, 1.9, 2.0, 10.0])
        
        # 2d interpolation
        self.f_cd0 = interpolate.interp1d(Mach_array,Cd0_array,kind='linear')
        
        return None
        
    def state2vecs_quat(self,u):
        # convert state vector u to vectors
        x = u[0:3]     # translation         :in fixed coordinate
        v = u[3:6]     # velocity            :in body coordinate
        q = u[6:10]    # atitude quaternion  :convert from fixed to body
        omega = u[10:] # angular velocity    :in body coordinate
    
        # convert to quaternion
        q2 = quaternion.as_quat_array(q)
    
        return x,v,q2,omega
        
                    
    def add_backup(self,t,u):
        # ==============================================
        # backup time, flag, state u
        # append vectors to the array "history"
        # ==============================================
               
        tmp = np.r_[t,self.flag,u]
        #try:
        self.backup = np.append(self.backup,[tmp],axis=0)
        #except:
        #    pass
        
        
    