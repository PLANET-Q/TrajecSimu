#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Sat Jul  8 23:06:17 2017

@author: shugo
"""
import sys
import numpy as np
from scipy.integrate import ode, odeint
import quaternion 

class simu_main():
    
    """
    ----------------------------------------------------
    ------     Methods for initial setup        --------
    ----------------------------------------------------
    """
    
    def __init__(self,params):
        # =============================================
        # this method is called when instance is created
        #
        # INPUT: params = contains all parameters needed
        # =============================================
        
        # assign parameters 
        self.set_all_param(params)
        
        
    def set_all_param(self,params):
        # =============================================
        # TO DO: assign parameters to the instance
        #
        # INPUT: params = contains all parameters needed (dict of dict)
        # =============================================
        
        # executive control
        try:
            self.dt = params['exec']['dt']              # time step
            self.t_max = params['exec']['t_max']        # maximum time
            self.N_record = params['exec']['N_record']  # record history every 20 iteration
            self.integ = params['exec']['integ']        # time integration scheme
        except:
            # display error message
            print 'Error: executive control info is missing'
            sys.exit()
        
        # launch condition
        try:
            # launcher property
            rail_length = params['launch']['rail_length']             # length of launch rail 
            self.elev_angle = params['launch']['elev_angle']          # angle of elevation [deg]
            self.azimuth = params['launch']['azimuth']                # north=0, east=90, south=180, west=270 [deg]
            self.rail_height = rail_length * np.sin(self.elev_angle * np.pi/180.) # height of launch rail in fixed coord.
            
            # air property
            self.T0 = params['launch']['T0']   # temperature [K] at 10m alt.
            self.p0 = params['launch']['p0']   # pressure [Pa] at 10m alt.
            
            # wind property
            wind_direction = params['launch']['wind_direction']  # azimuth where wind is blowing from [deg]
            angle_wind = (-wind_direction + 90.) * np.pi/180.    # angle to which wind goes (x orients east, y orients north)
            self.wind_unitvec = -np.array([np.cos(angle_wind), np.sin(angle_wind) ,0.])
            self.wind_speed = params['launch']['wind_speed']          # speed of wind [m/s] at 10m alt.
            self.Cwind = params['launch']['Cwind']                    # wind coefficient
        except:
            # display error message
            print 'Error: launch condition info is missing'
            sys.exit()
         
        # mass/inertia properties
        try:
            self.m_dry = params['rocket']['m_dry']        # dry weight of rocket i.e. exclude fuel
            self.m_fuel = params['rocket']['m_fuel']      # fule weight // NOTE: total weight at lift off = m_dry + m_fuel
            self.CG_dry = params['rocket']['CG_dry']      # CG of dried body (nose tip = 0)
            self.CG_fuel = params['rocket']['CG_fuel']    # CG of fuel (assume CG_fuel = const.)
            self.MOI_dry = params['rocket']['MOI_dry']    # dry MOI (moment of inertia)
            self.MOI_fuel = params['rocket']['MOI_fuel']  # dry MOI (moment of inertia)
        except:
            # display error message
            print 'Error: mass property info is missing'
            sys.exit()
            
        # aerodynamic properties
        try:
            self.CP_body = params['rocket']['CP_body']   # CP location without fins (budy tube+nose) (nose tip = 0)
            self.Cd0 = params['rocket']['Cd0']           # total 0 angle-of-attack drag coefficient
            self.X_area = params['rocket']['X_area']   # cross-sectional area
        
            # fin property
            self.LE_fins = params['rocket']['LE_fins']   # root leading edge of fin location (nose tip = 0)
            self.xFP = params['rocket']['xFP']  # array to store forcing-point locations (0=leading edge*root)
            self.fin_len = params['rocket']['fin_len']  # fin length
            self.r_arm = params['rocket']['r_arm'] # array of arm length for roll-spin
            self.alpha_fins = params['rocket']['alpha_fins']  # fin attachment angle
            self.dy_fin = params['rocket']['dy_fin']  # y step of fin discretization
        except:
            # display error message
            print 'Error: aerodynamic property info is missing'
            sys.exit()
        
        # thrust property
        try:
            self.t_MECO = params['engine']['t_MECO']
            self.thrustforce = params['engine']['thrust']
        except:
            # display error message
            print 'Error: engine property info is missing'
            
        # parachute property
        try:
            self.t_deploy = params['parachute']['t_deploy']
        except:
            # display error message
            print 'Error: parachute property info is missing'
            
    """
    ----------------------------------------------------
    ----     Method for main ODE integration      ------
    ----------------------------------------------------
    """
    
    def ODE_main(self):
        # =============================================
        # run ODE integration of rocket trajectory
        #
        # TO DO: ODE integration given all parameters
        #     
        # OUTPUT: result 
        # =============================================
        
        # ---------------------------------------------
        #      Initialization
        # ---------------------------------------------
        # set initial value
        t0 = 0.  # start time
        u0 = self.SetupICs() # initiate state vector

        # number of iteration
        self.N_iter = 0
            
        # initialize list to store the log
        #  histry = [times,flags,]
        self.history = [np.r_[t0,0,u0]]
    
        # set flag = 1 (on launch rail)
        self.flag = 1
    
        print '----------------------------'
        print '  We are GO for launch'
        print '----------------------------'
        print ' '
        
        
        if self.integ == 'lsoda_odeint':
            # ---------------------------------------------
            #      scipy.odeint
            # ---------------------------------------------
            # time array
            t = np.arange(t0,80,self.dt)
            try:
                # ---  run trajectory computation   ---
                sol = odeint(self.f_main, u0, t)
            except:
                pass
            
        else:
            # ---------------------------------------------
            #      scipy.ode
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
        #END IF
                
        print('----------------------------')
        print('   Completed trajectory computation ')
        print('----------------------------')
        print(' ')
        print(' ')
         
        
        
    """
    ----------------------------------------------------
    ----     Other methods for simulation         ------
    ----------------------------------------------------
    """

    
    def SetupICs(self):
        # =======================================
        # setup initial condition
        # return initial state vector u0: 13*1 array
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
        # right hand side of ODE
        #
        # INPUT: t    = time (scalar) 
        #        u    = state vector u (13*1)
        #
        # NOTE:  self.flag = 0 if before ignition
        #                    1 if on launch rail
        #                    2 if thrusted flight
        #                    3 if inertial flight
        # =======================================
        
        # swap input when we use scipy.ubtegrate
        if self.integ != 'lsoda_odeint':
            tmp = u
            u = t
            t = tmp
        #END IF
    
        # count iteration
        self.N_iter += 1
        
        # record history
        if np.mod(self.N_iter, self.N_record) == 0:
            self.record_history(t,u)
        #END IF
        
        # --------------------------
        #   subtract vectors
        # --------------------------
        # x =     translation         :in fixed coordinate
        # v =     velocity            :in body coordinate
        # q =     atitude quaternion  :convert from fixed to body
        # omega = angular velocity    :in body coordinate
        x,v,q,omega = self.state2vecs_quat(u)  
        
        
        # ----------------------------
        #   flight mode
        # ----------------------------
        if self.flag==1 and x[2] > self.rail_height:
            # detect launch clear
            print('----------------------------')
            print('  Launcher-clear at t = ',t,'[s]')
            print('  current speed: ',np.linalg.norm(v),'[m/s]')
            print('----------------------------')
            print(' ')
            
            # record history
            self.record_history(t,u)
            # switch into thrusted flight
            self.flag = 2   
            
        elif (self.flag==1 or self.flag==2) and t >= self.t_MECO:
            # detect MECO
            print('----------------------------')
            print('  MECO at t = ',t,'[s]')
            print('  current altitude: ',x[2],'[m]')
            print('          speed:    ',np.linalg.norm(v),'[m/s]')
            print('----------------------------')
            print(' ')
            
            # record history
            self.record_history(t,u)
            # switch into inertial flight
            self.flag = 3   
            
        elif (self.flag==2 or self.flag==3) and t >= self.t_deploy:
            # detect parachute deployment
            print('----------------------------')
            print('  Parachute deployed at t = ',t,'[s]')
            print('  current altitude: ',x[2],'[m]')
            print('          speed:    ',np.linalg.norm(v),'[m/s]')
            print('----------------------------')
            print(' ')
            
            # record history
            self.record_history(t,u)
            # switch into parachute fall
            self.flag = 4   
            
        elif self.flag > 1 and x[2] < 0.:
            # detect landing
            print('----------------------------')
            print('  Landing at t = ',t,'[s]')
            print('  landing speed: ',np.linalg.norm(v),'[m/s]')
            print('          location x = ',x[0],'[m]')
            print('                   y = ',x[1],'[m]')
            print('----------------------------')
            print(' ')
            
            # record history
            self.record_history(t,u)
            # landed
            self.flag = 5 
            # quit integration
            if self.integ == 'lsoda_odeint':
                return 0
                
            
        #END IF 
            
            
        # call mass properties
        mass,MOI,d_dt_MOI,CG = self.mass_MOI(t)
            
        # ----------------------------
        #    Direction Cosine Matrix for input q
        # ----------------------------
        # Tbl = transform from local(fixed) coordinate to body coord.
        #     note: as_rotation_matrix is for vector rotation
        #         -> for coordinate rotation, input conj(q)
        Tbl= quaternion.as_rotation_matrix(np.conj(q))
        
        
        # ----------------------------
        #    1. Translation
        # ----------------------------
        # translation time rate = velocity
        # convert v from body coord. to fixed coord.
        dx_dt = np.dot(Tbl.T,v)
        
        
        # ----------------------------
        #    2. Velocity
        # ----------------------------
        # velocity time rate = mass * force
        # force = rotation effect + aero + thrust + gravity 
        #      where  aero   = aero(u):   function of state u
        #             thrust = thrust(t): function of time t
        #             mass   = mass(t):   function of time t   
    
        # call aerodynamic force/moment
        aeroF,aeroM = self.aero(x,v,omega,Tbl,CG)     
        
        grav = np.array([0.,0.,-9.81])  # gravitational accel. in fixed coord.
        
        # set external force depending on the state
        if self.flag == 1:
            # on launch rail -> du/dx only. Consider rail-rocket friction 
            # total acceleration 
            dv_dt = -np.cross(omega,v) + np.dot(Tbl,grav) + (aeroF + self.thrust(t) + self.friction() ) / mass
            # cancell out y,z
            dv_dt = np.array([dv_dt[0],0.,0.])
            
        elif self.flag == 2:
            # thrust ON
            # total acceleration 
            dv_dt = -np.cross(omega,v) + np.dot(Tbl,grav) + (aeroF + self.thrust(t)) / mass
            
        elif self.flag == 3 or self.flag == 5:
            # thrust OFF
            # total acceleration 
            dv_dt = -np.cross(omega,v) + np.dot(Tbl,grav) + aeroF / mass
        #END IF
        
        # ----------------------------
        #    3. Atitude 
        # ----------------------------
        # represented in quaternion which rotates from fixed to body
        # quaternion differential equation
        
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
        
        if self.flag == 1:
            # on launch rail -> no angular velocity change
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
        #    Set back in state vector form
        # ----------------------------
        du_dt = np.r_[dx_dt,dv_dt,dq_dt,domega_dt]
        
        return du_dt
        
    def mass_MOI(self,t):
        # =======================================
        # returns mass properties of rocket
        # 
        # INPUT:  t = time
        # OUTPUT: mass = total mass of rocket
        #         MOI = total moment of inertia wrt CG
        #         d_dt_MOI = MOI time rate
        #         CG = center of gravity location from the nose tip
        # =======================================

        if t >= self.t_MECO:
            # mass
            mass = self.m_dry
            # moment of inertia
            MOI = self.MOI_dry
            # time rate of MOI
            d_dt_MOI = np.zeros(3)
            # total CG location
            CG = self.CG_dry
            
        else:
            # mass
            mass = self.m_dry + (self.t_MECO-t)/self.t_MECO * self.m_fuel
    
            # fuel comsumption rate (assumed linear)
            r = (self.t_MECO-t)/self.t_MECO 
            # total CG location
            CG = (self.CG_dry*self.m_dry + self.CG_fuel*r*self.m_fuel) / (self.m_dry + r*self.m_fuel)
            # total MOI using parallel axis theorem
            tmp = np.array([0.,1.,1.])
            MOI = self.MOI_dry + tmp*self.m_dry*(CG-self.CG_dry)**2. + r*self.MOI_fuel + tmp*(CG-self.CG_fuel)*(r*self.m_fuel)**2.
    
            # finite differencing
            h = 1.E-5
            r2 = (self.t_MECO-t+h)/self.t_MECO 
            # total CG location
            CG2 = (self.CG_dry*self.m_dry + self.CG_fuel*r2*self.m_fuel) / (self.m_dry + r2*self.m_fuel)
            # total MOI using parallel axis theorem
            tmp = np.array([0.,1.,1.])
            MOI2 = self.MOI_dry + tmp*self.m_dry*(CG2-self.CG_dry)**2. + r2*self.MOI_fuel + tmp*(CG2-self.CG_fuel)*(r2*self.m_fuel)**2.
    
            # finite differencing
            r3 = (self.t_MECO-t-h)/self.t_MECO 
            # total CG location
            CG3 = (self.CG_dry*self.m_dry + self.CG_fuel*r3*self.m_fuel) / (self.m_dry + r3*self.m_fuel)
            # total MOI using parallel axis theorem
            tmp = np.array([0.,1.,1.])
            MOI3 = self.MOI_dry + tmp*self.m_dry*(CG3-self.CG_dry)**2. + r3*self.MOI_fuel + tmp*(CG3-self.CG_fuel)*(r3*self.m_fuel)**2.
    
            d_dt_MOI = (MOI2 - MOI3) / (2*h)
        #END IF
        
        return mass, MOI, d_dt_MOI, CG
                
    
    def thrust(self,t):
        # =======================================
        # returns thrust force 
        # 
        # INPUT:  t = time
        # OUTPUT: T = thrust vector in body coord.
        # =======================================
        T = np.array([self.thrustforce,0.,0.])
        return T
        
        
    def friction(self):
        # =======================================
        # returns friction force (rail-rocket contact) 
        # 
        # INPUT:  t = time
        # OUTPUT: fric = friction force vector in body coord.
        # =======================================
        return np.array([-5.,0.,0.])
            
        
            
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
        # air velocity = -rocket_velocity + wind_velocity
        #    NOTE: wind(x) is wind velocity in local coord. need conversion to body coord.
        v_air = -v + np.dot( Tbl,self.wind(x[2]) )   
        air_speed = np.linalg.norm(v_air) # air speed (positive scalar)
    
        # total angle-of-attack
        alpha = np.arccos( -v_air[0]/air_speed )
        if np.isnan(alpha):
            alpha = 0.
        
        # roll-angle
        phi = np.arctan( v_air[1]/v_air[2] )
        if np.isnan(phi):
            phi = np.arctan( v_air[1]+0.0001/v_air[2]+0.0001 )
            
        
        # ----------------------------------
        #   force on body excluding fins
        # ----------------------------------
        # air property at the altitude
        rho,a = self.standard_air(x[2])
        
        # Mach number
        Mach = air_speed/a
    
        # drag/lift coefficient of body
        Cd,Cl = self.rocket_coeff_nofins(Mach,alpha)
        
        # convert coefficient to body coord.
        cosa = np.cos(alpha)
        sina = np.sin(alpha)
        #C1b = -Cl*sina + Cd*cosa
        #C2b = -(Cl*cosa + Cd*sina)*np.sin(phi)
        #C3b = (Cl*cosa + Cd*sina)*np.cos(phi)
        C = np.array([- (-Cl*sina + Cd*cosa), \
                      -(Cl*cosa + Cd*sina)*np.sin(phi), \
                      -(Cl*cosa + Cd*sina)*np.cos(phi)])    # need check here
        
        # force on CP of body
        force_body = 0.5 * rho * air_speed**2. * self.X_area * C
        
        # moment generated by body wrt CG
        moment_body = np.cross( np.array([CG-self.CP_body,0.,0.]) , force_body )
        
        # ----------------------------------
        #   force on fins
        # ----------------------------------
        # force and moment (wrt fin Leading edge)
        force_fin, moment_fin_tmp = self.fin_aero(v_air,omega[0],rho)
        # moment generated by fin wrt CG
        moment_fin = moment_fin_tmp + np.cross( np.array([CG-self.LE_fins,0.,0.]) , force_fin )
        
        # ----------------------------------
        #  total force and moment
        # ----------------------------------
        # total force
        force_all = force_body + force_fin
        # total moment
        moment_all = moment_body + moment_fin
        
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
    
        
        
    def standard_air(self,h):
        # ==============================================
        # returns air property given an altitude 
        # INPUT: h = altitude [m]
        # ==============================================
        
        # temperature goes down 0.0065K/m until it reaches -56.5C (216.5K)
        #                                       it is approximately 11km
        T = self.T0 - 0.0065*h # [K]
        if T < 216.5:
            # temperature is const at 216.5 K for alt. < 20km
            T = 216.5
            
        # pressure
        p = self.p0 * (T/self.T0)**5.256  #[Pa]

        # density
        rho = p/(287.15*T) #[kg/m^3]
        
        # acoustic speed
        a = np.sqrt(1.4*287.15*T) # [m/s]
        
        return rho,a 
        
    
    def wind(self,h):
        # ==============================================
        # returns wind vector given an altitude
        #  follows "wind profile power law"
        #
        # INPUT: h = altitude [m]
        # ==============================================
        
        # wind velocity in 
        wind_vec = self.wind_unitvec * self.wind_speed * (h/10.)**self.Cwind  

        return wind_vec
        
        
        
        
    def rocket_coeff_nofins(self, Mach,alpha):
        
        # lift coefficient slope
        k = 1.0
        
        # drag coefficient slope
        k2 = 1.0
        
        if Mach > 0.8:
            Mach = 0.8
        
        Cl = k * alpha / np.sqrt(1-Mach**2.)
    
        Cd = self.Cd0 + k2*alpha**2. / np.sqrt(1-Mach**2.)
    
        return Cd, Cl
        
        
        
        
    def state2vecs_quat(self,u):
        # convert state vector u to vectors
        x = u[0:3]     # translation         :in fixed coordinate
        v = u[3:6]     # velocity            :in body coordinate
        q = u[6:10]    # atitude quaternion  :convert from fixed to body
        omega = u[10:] # angular velocity    :in body coordinate
    
        # convert to quaternion
        q2 = quaternion.as_quat_array(q)
    
        return x,v,q2,omega
        
        
                    
    def record_history(self,t,u):
        # ==============================================
        # record history
        # append vectors to the array "history"
        # ==============================================
               
        tmp = np.r_[t,self.flag,u]
        #try:
        self.history = np.append(self.history,[tmp],axis=0)
        #except:
        #    pass
        
        
    