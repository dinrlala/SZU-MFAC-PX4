/****************************************************************************
 *
 *   Copyright (c) 2018 - 2019 PX4 Development Team. All rights reserved.
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions
 * are met:
 *
 * 1. Redistributions of source code must retain the above copyright
 *    notice, this list of conditions and the following disclaimer.
 * 2. Redistributions in binary form must reproduce the above copyright
 *    notice, this list of conditions and the following disclaimer in
 *    the documentation and/or other materials provided with the
 *    distribution.
 * 3. Neither the name PX4 nor the names of its contributors may be
 *    used to endorse or promote products derived from this software
 *    without specific prior written permission.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
 * "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
 * LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS
 * FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE
 * COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT,
 * INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING,
 * BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS
 * OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED
 * AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
 * LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN
 * ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
 * POSSIBILITY OF SUCH DAMAGE.
 *
 ****************************************************************************/

/**
 * @file PositionControl.cpp
 */

#include "PositionControl.hpp"
#include "ControlMath.hpp"
#include <float.h>
#include <mathlib/mathlib.h>
#include <matrix/matrix/math.hpp>
#include <px4_platform_common/defines.h>
#include <geo/geo.h>
#include <px4_platform_common/log.h>
//uORB
#include <uORB/uORB.h>
#include <uORB/topics/trajectory_setpoint.h>
#include <uORB/topics/vehicle_attitude_setpoint.h>
#include <uORB/topics/vehicle_local_position_setpoint.h>
#include <uORB/topics/vehicle_status.h>
#include <uORB/topics/vehicle_land_detected.h>

//#include <Eigen/Core>

using namespace matrix;

//enum class MFACState { PID_INIT, MFAC_ACTIVE };
//static MFACState _mfac_state = MFACState::PID_INIT;
const float eps = 1e-6f;
static int print_counter = 0;
//static int update_counter = 5;

static inline void mfac_shift_histories(
	matrix::Matrix<float,3,3> &uk,
	matrix::Matrix<float,3,3> &ek)
{
	for (int i = 0; i < 3; i++) {
		uk(2,i) = uk(1,i);
		uk(1,i) = uk(0,i);
		ek(2,i) = ek(1,i);
		ek(1,i) = ek(0,i);
	}
}

static inline void mfac_init_histories(
	matrix::Matrix<float,3,3> &uk,
	matrix::Matrix<float,3,3> &ek)
{
	for (int i = 0; i < 3; i++) {
		uk(2,i) = uk(1,i) = uk(0,i);
		ek(2,i) = ek(1,i) = ek(0,i);
	}
}


//提取某一列的thetack和thetack1
Vector3f get_thetack(const Matrix<float, 6, 3>& mat, size_t col_index) {
    	matrix::Vector3f result;
    	for (size_t i = 0; i < 3; ++i) {
        	result(i) = mat(i, col_index);
    	}
    	return result;
	}

Vector3f get_thetack1(const Matrix<float, 6, 3>& mat, size_t col_index) {
    	matrix::Vector3f result;
    	for (size_t i = 0; i < 3; ++i) {
        	result(i) = mat(3 + i, col_index);  // 从第3行开始
    	}
    	return result;
}

void edit_thetack(Matrix<float, 6, 3>& targetMatrix, size_t col_index, const Vector3f& targetVector) {
    for (size_t i = 0; i < 3; ++i) {
        targetMatrix(i, col_index) = targetVector(i);
    }
}



const trajectory_setpoint_s PositionControl::empty_trajectory_setpoint = {0, {NAN, NAN, NAN}, {NAN, NAN, NAN}, {NAN, NAN, NAN}, {NAN, NAN, NAN}, NAN, NAN};

void PositionControl::setVelocityGains(const Vector3f &P, const Vector3f &I, const Vector3f &D)
{
	_gain_vel_p = P;
	_gain_vel_i = I;
	_gain_vel_d = D;
}

void PositionControl::setVelocityGainsMFAC(const Vector3f &LAMBDAC, const Vector3f &LAMBDAM, const Vector3f &THETACTEMPXY,const Vector3f &THETACTEMPZ,const Vector3f &THETAM)//设置MFAC要用到的参数，在参数表中直接读取
{
	_mfac_vel_lambdac = LAMBDAC;
	_mfac_vel_lambdam = LAMBDAM;
	_mfac_vel_thetac_temp_xy = THETACTEMPXY;
	_mfac_vel_thetac_temp_z = THETACTEMPZ;
	_mfac_vel_thetam = THETAM;
}

void PositionControl::setControllerMode(const int &MODE)//设置使用什么Controller
{
	_vel_con_choose = MODE;
}

void PositionControl::setMFACMode(const int &MODE)//设置是否固定增益
{
	_mfac_sol_mode = MODE;
}

void PositionControl::setMFACPIDInit(const int &TIME)//设置PID初始化时间
{
	_mfac_pid_init = TIME;
}

void PositionControl::setXYThetac2Limit(const float &LIMIT)//设置thetac限制
{
	_mfac_vel_thetac_xy_thetac2Limit = LIMIT;
}


void PositionControl::setXYThetac3Limit(const float &LIMIT3)//设置thetac限制
{
	_mfac_vel_thetac_xy_thetac3Limit = LIMIT3;
}


void PositionControl::setXYAccLimit(const float &LIMITACC)//设XY加速度限制
{
	_mfac_vel_acc_xy_Limit = LIMITACC;
}





void PositionControl::setVelocityLimits(const float vel_horizontal, const float vel_up, const float vel_down)
{
	_lim_vel_horizontal = vel_horizontal;
	_lim_vel_up = vel_up;
	_lim_vel_down = vel_down;
}

void PositionControl::setThrustLimits(const float min, const float max)
{
	// make sure there's always enough thrust vector length to infer the attitude
	_lim_thr_min = math::max(min, 10e-4f);
	_lim_thr_max = max;
}

void PositionControl::setHorizontalThrustMargin(const float margin)
{
	_lim_thr_xy_margin = margin;
}

void PositionControl::updateHoverThrust(const float hover_thrust_new)
{
	// Given that the equation for thrust is T = a_sp * Th / g - Th
	// with a_sp = desired acceleration, Th = hover thrust and g = gravity constant,
	// we want to find the acceleration that needs to be added to the integrator in order obtain
	// the same thrust after replacing the current hover thrust by the new one.
	// T' = T => a_sp' * Th' / g - Th' = a_sp * Th / g - Th
	// so a_sp' = (a_sp - g) * Th / Th' + g
	// we can then add a_sp' - a_sp to the current integrator to absorb the effect of changing Th by Th'
	const float previous_hover_thrust = _hover_thrust;
	setHoverThrust(hover_thrust_new);

	_vel_int(2) += (_acc_sp(2) - CONSTANTS_ONE_G) * previous_hover_thrust / _hover_thrust
		       + CONSTANTS_ONE_G - _acc_sp(2);
}

void PositionControl::setState(const PositionControlStates &states)
{
	_pos = states.position;
	_vel = states.velocity;
	_yaw = states.yaw;
	_vel_dot = states.acceleration;
}

void PositionControl::setInputSetpoint(const trajectory_setpoint_s &setpoint)
{
	_pos_sp = Vector3f(setpoint.position);
	_vel_sp = Vector3f(setpoint.velocity);
	_acc_sp = Vector3f(setpoint.acceleration);
	_yaw_sp = setpoint.yaw;
	_yawspeed_sp = setpoint.yawspeed;
}

bool PositionControl::update(const float dt)
{
	bool valid = _inputValid();

	if (valid) {
		_positionControl();
		switch (_vel_con_choose)
		{
		case 0: //位置PID
			_velocityControl(dt);
			break;
		case 1: //CDL-MFAC
			_velocityControlMFAC(dt);
			break;
		default://默认PID
			_velocityControl(dt);
			break;
		}

		_yawspeed_sp = PX4_ISFINITE(_yawspeed_sp) ? _yawspeed_sp : 0.f;
		_yaw_sp = PX4_ISFINITE(_yaw_sp) ? _yaw_sp : _yaw; // TODO: better way to disable yaw control
	}

	// There has to be a valid output acceleration and thrust setpoint otherwise something went wrong
	return valid && _acc_sp.isAllFinite() && _thr_sp.isAllFinite();
}

void PositionControl::_positionControl()
{
	// P-position controller
	Vector3f vel_sp_position = (_pos_sp - _pos).emult(_gain_pos_p);
	// Position and feed-forward velocity setpoints or position states being NAN results in them not having an influence
	ControlMath::addIfNotNanVector3f(_vel_sp, vel_sp_position);
	// make sure there are no NAN elements for further reference while constraining
	ControlMath::setZeroIfNanVector3f(vel_sp_position);

	// Constrain horizontal velocity by prioritizing the velocity component along the
	// the desired position setpoint over the feed-forward term.
	_vel_sp.xy() = ControlMath::constrainXY(vel_sp_position.xy(), (_vel_sp - vel_sp_position).xy(), _lim_vel_horizontal);
	// Constrain velocity in z-direction.
	_vel_sp(2) = math::constrain(_vel_sp(2), -_lim_vel_up, _lim_vel_down);
}

void PositionControl::_velocityControl(const float dt)
{
	// Constrain vertical velocity integral
	_vel_int(2) = math::constrain(_vel_int(2), -CONSTANTS_ONE_G, CONSTANTS_ONE_G);

	// PID velocity control
	Vector3f vel_error = _vel_sp - _vel;
	Vector3f acc_sp_velocity = vel_error.emult(_gain_vel_p) + _vel_int - _vel_dot.emult(_gain_vel_d);//这里分了PID三段



	// No control input from setpoints or corresponding states which are NAN
	ControlMath::addIfNotNanVector3f(_acc_sp, acc_sp_velocity);

	_accelerationControl();

	// Integrator anti-windup in vertical direction
	if ((_thr_sp(2) >= -_lim_thr_min && vel_error(2) >= 0.f) ||
	    (_thr_sp(2) <= -_lim_thr_max && vel_error(2) <= 0.f)) {
		vel_error(2) = 0.f;
	}

	// Prioritize vertical control while keeping a horizontal margin
	const Vector2f thrust_sp_xy(_thr_sp);
	const float thrust_sp_xy_norm = thrust_sp_xy.norm();
	const float thrust_max_squared = math::sq(_lim_thr_max);

	// Determine how much vertical thrust is left keeping horizontal margin
	const float allocated_horizontal_thrust = math::min(thrust_sp_xy_norm, _lim_thr_xy_margin);
	const float thrust_z_max_squared = thrust_max_squared - math::sq(allocated_horizontal_thrust);

	// Saturate maximal vertical thrust
	_thr_sp(2) = math::max(_thr_sp(2), -sqrtf(thrust_z_max_squared));

	// Determine how much horizontal thrust is left after prioritizing vertical control
	const float thrust_max_xy_squared = thrust_max_squared - math::sq(_thr_sp(2));
	float thrust_max_xy = 0.f;

	if (thrust_max_xy_squared > 0.f) {
		thrust_max_xy = sqrtf(thrust_max_xy_squared);
	}

	// Saturate thrust in horizontal direction
	if (thrust_sp_xy_norm > thrust_max_xy) {
		_thr_sp.xy() = thrust_sp_xy / thrust_sp_xy_norm * thrust_max_xy;
	}

	// Use tracking Anti-Windup for horizontal direction: during saturation, the integrator is used to unsaturate the output
	// see Anti-Reset Windup for PID controllers, L.Rundqwist, 1990
	const Vector2f acc_sp_xy_produced = Vector2f(_thr_sp) * (CONSTANTS_ONE_G / _hover_thrust);
	const float arw_gain = 2.f / _gain_vel_p(0);

	// The produced acceleration can be greater or smaller than the desired acceleration due to the saturations and the actual vertical thrust (computed independently).
	// The ARW loop needs to run if the signal is saturated only.
	const Vector2f acc_sp_xy = _acc_sp.xy();
	const Vector2f acc_limited_xy = (acc_sp_xy.norm_squared() > acc_sp_xy_produced.norm_squared())
					? acc_sp_xy_produced
					: acc_sp_xy;
	vel_error.xy() = Vector2f(vel_error) - arw_gain * (acc_sp_xy - acc_limited_xy);

	// Make sure integral doesn't get NAN
	ControlMath::setZeroIfNanVector3f(vel_error);
	// Update integral part of velocity control
	_vel_int += vel_error.emult(_gain_vel_i) * dt;
}



//速度环MFAC控制器
void PositionControl::_velocityControlMFAC(const float dt)// dt为时间步长
{
	//const bool armed = (_vehicle_status.arming_state == vehicle_status_s::ARMING_STATE_ARMED);
	//const bool landed = _land_detected.landed;
	Vector3f vel_error = _vel_sp - _vel;// 期望速度-实际速度
	const float vel_err_thresh = 0.05f;   // m/s
	const float vel_sp_thresh  = 0.05f;

	if (!_mfac_allow) {
        	//未解锁不进入初始化阶段
        	_velocityControl(dt);
        	return;
    	}


	//if (landed) {
   	// 	_mfac_state = MFACState::PID_INIT;

    	//	Vector3f acc_sp_velocity =
        //		vel_error.emult(_gain_vel_p)
        //		+ _vel_int
        //		- _vel_dot.emult(_gain_vel_d);

    	//	ControlMath::addIfNotNanVector3f(_acc_sp, acc_sp_velocity);
    	//return;
	//}

	//入口保护，防止非法的速度指令
	if ((!_vel.isAllFinite() || !_vel_sp.isAllFinite()) && _mfac_state == MFACState::PID_INIT) {
    		return;
	}


	//还是防止非法的速度指令
	if (!_velk1.isAllFinite()) {
    		_velk1 = _vel;
	}
	//保护性参数，配合限幅保护用
	//const float thetac0_max = 20.0f;         // thetak0最大绝对值，经验值
	//const float thetac1_max = _mfac_vel_thetac_xy_thetac2Limit > 0 ?
        //_mfac_vel_thetac_xy_thetac2Limit : 5.0f; // fallback
	//const float thetac2_max = _mfac_vel_thetac_xy_thetac3Limit > 0 ?
        //_mfac_vel_thetac_xy_thetac3Limit : 5.0f;

	if(_vel_con_choose==1){//省空间用
	//这里thetack定义为一个6x3的矩阵，给未来的自己和后人，勿看也看不懂
	//   x轴               y轴                  z轴
	//  thetak_x(0)    thetak_y(0)        thetak_z(0)
	//  thetak_x(1)    thetak_y(1)        thetak_z(1)
	//  thetak_x(2)    thetak_y(2)        thetak_z(2)
	//  thetak1_x(0)    thetak1_y(0)        thetak1_z(0)
	//  thetak1_x(1)    thetak1_y(1)        thetak1_z(1)
	//  thetak1_x(2)    thetak1_y(2)        thetak1_z(2)

	//速度误差历史更新

	for(int i=0;i<3;i++){
		ek(0,i)=_vel_sp(i)-_vel(i);
	}
	//thetack拼接成3X3矩阵
	//for (int i = 0; i < 2; i++) {
    	//	_mfac_vel_thetac(0,i) = _mfac_vel_thetac_temp_xy(0);
    	//	_mfac_vel_thetac(1,i) = _mfac_vel_thetac_temp_xy(1);
    	//	_mfac_vel_thetac(2,i) = _mfac_vel_thetac_temp_xy(2);
	//}
	//_mfac_vel_thetac(0,2) = _mfac_vel_thetac_temp_z(0);
    	//_mfac_vel_thetac(1,2) = _mfac_vel_thetac_temp_z(1);
    	//_mfac_vel_thetac(2,2) = _mfac_vel_thetac_temp_z(2);


	if(_mfac_state == MFACState::PID_INIT){//如果未初始化，使用PID进行初始化
		//初始化输入输出误差值，使用PID进行

		Vector3f acc_sp_velocity = vel_error.emult(_gain_vel_p) + _vel_int - _vel_dot.emult(_gain_vel_d);

		ControlMath::addIfNotNanVector3f(_acc_sp, acc_sp_velocity);
		count++;
		for (int i = 0; i < 3; i++) {
                    	uk(0, i) = acc_sp_velocity(i);
               	}
		mfac_init_histories(uk, ek);
		if(count>=_mfac_pid_init){//初始化完毕
                	//初始化thetam
                	for (int i = 0; i < 3; i++) {
                    		thetamk(0, i) = _mfac_vel_thetam(i);
                    		thetamk(1, i) = _mfac_vel_thetam(i);
               		 }

                	// 把参数填到 thetack 的两份历史（rows 0..2 和 3..5）
                	//初始化thetac
                	for (int j = 0; j < 3; j++) {
                    		thetack(j + 0, 0) = _mfac_vel_thetac_temp_xy(j);
                    		thetack(j + 3, 0) = _mfac_vel_thetac_temp_xy(j);
				thetack_init(j + 0, 0) = _mfac_vel_thetac_temp_xy(j);
                    		thetack_init(j + 3, 0) = _mfac_vel_thetac_temp_xy(j);
               	 	}
                	for (int j = 0; j < 3; j++) {
                    		thetack(j + 0, 1) = _mfac_vel_thetac_temp_xy(j);
                    		thetack(j + 3, 1) = _mfac_vel_thetac_temp_xy(j);
				thetack_init(j + 0, 1) = _mfac_vel_thetac_temp_xy(j);
                    		thetack_init(j + 3, 1) = _mfac_vel_thetac_temp_xy(j);
                	}
                	for (int j = 0; j < 3; j++) {
                    		thetack(j + 0, 2) = _mfac_vel_thetac_temp_z(j);
                    		thetack(j + 3, 2) = _mfac_vel_thetac_temp_z(j);
				thetack_init(j + 0, 2) = _mfac_vel_thetac_temp_z(j);
                    		thetack_init(j + 3, 2) = _mfac_vel_thetac_temp_z(j);
               		 }
                	//uk 的最新值为当前 PID 输出
                _mfac_state = MFACState::MFAC_ACTIVE;
                count = 0;
		}
	}

	if(_mfac_state == MFACState::MFAC_ACTIVE){
		//thetac和thetam的历史更新
		for(int i = 0;i < 3;i++){
			thetamk(1,i)=thetamk(0,i);
			for(int j = 0;j<3;j++){
				thetack(j+3,i)=thetack(j,i);
			}
	 	}
		//更新历史信息ek,uk
		mfac_shift_histories(uk, ek);

		if(!_mfac_sol_mode){//自适应模式
			//求Hk
			for (int j = 0; j < 3; j++) {
        			Hk(0,j) = -ek(0,j);
    			}
    			for (int j = 0; j < 3; j++) {
        			Hk(1,j) = ek(0,j) - ek(1,j);
    			}
    			for (int j = 0; j < 3; j++) {
        			Hk(2,j) = ek(1,j) - ek(2,j);
    			}
			//求XY的thetamk
			for(int i =0;i<2;i++){
				float thetamTemp = _mfac_vel_lambdam(i)+(uk(1,i)-uk(2,i))*(uk(1,i)-uk(2,i));
				thetamk(0,i)=thetamk(1,i)+(_vel(i)-_velk1(i)-(thetamk(1,i)*(uk(1,i)-uk(2,i))))*(uk(1,i)-uk(2,i))/thetamTemp;
			}
			//求XY的thetack
			for(int i=0;i<2;i++){
				matrix::Vector3f thetacik,thetacik1;
				thetacik1 = get_thetack1(thetack,i);
				matrix::Vector3f Hk_col = matrix::Vector3f(Hk.col(i));
				float Hknorm = Hk_col.norm();
				float tempThetac=(get_thetack1(thetack,i).transpose()*Hk_col)(0, 0);
				thetacik = thetacik1 + (thetamk(0,i)*Hk_col*(_vel_sp(i)-_vel(i)-thetamk(0,i)*tempThetac))/(_mfac_vel_lambdac(i)+Hknorm*Hknorm);
				//**************保护性措施***************************************
				// 限幅保护
    				//if (fabsf(thetacik(0)) > thetac0_max) thetacik(0) = sign(thetacik(0)) * thetac0_max;
    				//if (fabsf(thetacik(1)) > thetac1_max) thetacik(1) = sign(thetacik(1)) * thetac1_max;
    				//if (fabsf(thetacik(2)) > thetac2_max) thetacik(2) = sign(thetacik(2)) * thetac2_max;
   				 //主元符号一致性保护
    				//float Hk0 = Hk_col(0);
    				//if (fabsf(Hk0) > eps) {
        				// 期望 thetac(0) * Hk0 > 0 (因为 uktemp = thetac^T * Hk，想让主项为正贡献)
        			//	if (thetacik(0) * Hk0 < 0.f) {
         			   	// 将主元强制为与 Hk0 同号的保守值（不直接取反大幅度变更）
          			// 	float safe_val = 0.05f * sign(Hk0); // 小值修正，防止突变
           			// 	thetacik(0) = safe_val;
       			 	//	}
				//}
   				//平滑更新保护
    				// 将新 thetacik 和 当前 thetac( rows 0..2 ) 做平滑融合，避免一次性跳变
    				const float alpha = 0.6f; // 新增量权重（0..1），越小更新越保守
    				matrix::Vector3f thetac_current = get_thetack(thetack, i); // rows 0..2
    				matrix::Vector3f thetac_new = thetac_current * (1.0f - alpha) + thetacik * alpha;
				edit_thetack(thetack,i,thetac_new);
			}
			// 推力异常处理
			// 在调用本段前确保 Hk 已经计算，uk/ek 已 shift，thetamk(0/1) 可用

			const float vel_delta_thresh = 0.01f;    // 速度变化很小时不要更新（单位 m/s）
			const float thrust_sat_ratio = 0.90f;    // 当垂直推力接近这个比例认为受限
			 // fallback
			const float thetac_forget = 0.950f;      // 遗忘因子（<1 会缓慢衰减旧的thetac，防止长期累积），贴近1

			// 计算当前垂直推力是否接近饱和，饱和则冻结更新
			const float thrust_z_norm = fabsf(_thr_sp(2));
			const float thrust_z_max = _lim_thr_max;
			bool thrust_z_near_sat = (thrust_z_max > 1e-6f) && (thrust_z_norm > thrust_sat_ratio * thrust_z_max);

			// 计算 Hk 范数（用于判断是否有激励）
			for (int col = 0; col < 2; col++) {
			//    matrix::Vector3f Hk_col = matrix::Vector3f(Hk.col(col));
			//    float Hk_norm = Hk_col.norm();
			//const float Hk_norm_thresh = 1e-4f;      // Hk 太小表示激励不足 -> 跳过更新
    			// 若 Hk 激励不足，跳过本列的自适应更新（保留历史 thetac）
   			// if (Hk_norm < Hk_norm_thresh) {
       			 // 轻微遗忘，避免长期累积（保持一定衰减）
        		//	for (int r = 0; r < 3; r++) {
        		//	    thetack(r, col) *= thetac_forget;
        		//	}
        		//	continue;
    			//}

    			// 若垂直推力接近饱和，禁止 XY 自适应（避免学到被削弱的映射）
    			if (thrust_z_near_sat) {
        			for (int r = 0; r < 3; r++) {
            			thetack(r, col) = thetac_forget*thetack(r, col)+(1-thetac_forget)*thetack_init(r,col); // 只做遗忘，不做新学习
        			}
       			 continue;
    			}
    			// 若速度几乎没变化，跳过更新（避免噪声/积分导致学错）
   			 if (fabsf(_vel(col) - _velk1(col)) < vel_delta_thresh) {
        			for (int r = 0; r < 3; r++) {
        			    thetack(r, col)=thetac_forget*thetack(r, col)+(1-thetac_forget)*thetack_init(r,col);
        			}
        			continue;
    			}
			}
			//误差变化量的限制,第二项
			//for(int i=0;i<2;i++){
			//	matrix::Vector3f thetacktemp=get_thetack(thetack,i);
			//	if(thetacktemp(1)>_mfac_vel_thetac_xy_thetac2Limit||thetacktemp(1)<-_mfac_vel_thetac_xy_thetac2Limit)
			//	thetacktemp(1)=sign(thetacktemp(1))*_mfac_vel_thetac_xy_thetac2Limit;
			//	edit_thetack(thetack,i,thetacktemp);
			//}
			//误差变化量的限制,第三项
			//for(int i=0;i<2;i++){
			//	matrix::Vector3f thetacktemp=get_thetack(thetack,i);
			//	if(thetacktemp(2)>_mfac_vel_thetac_xy_thetac3Limit||thetacktemp(2)<-_mfac_vel_thetac_xy_thetac3Limit)
			//	thetacktemp(2)=sign(thetacktemp(2))*_mfac_vel_thetac_xy_thetac3Limit;
			//	edit_thetack(thetack,i,thetacktemp);
			//}
			//求XY的uk
			//float acc_change_limit = 0.01f; // 每步允许的最大加速度变化量
			for(int i=0;i<2;i++){
				//float acc_delta = uk(0, i) - uk(1, i);
				matrix::Vector3f Hk_col = matrix::Vector3f(Hk.col(i));
				float uktemp = (get_thetack(thetack,i).transpose()*Hk_col)(0, 0);
				uk(0,i)=uk(1,i)+uktemp;
			//更新uk
				bool has_excitation =
    					fabsf(_vel_sp(i)) > vel_sp_thresh ||
    					fabsf(_vel(i))    > vel_err_thresh;
				if (!has_excitation) {
    				//判断为自激励则冻结thetack和u,只允许遗忘。
    					uk(0,i) = uk(1,i);
    					for (int row = 0; row < 6; row++) {
   					 	thetack(row, i) = thetac_forget * thetack(row, i) + (1.f - thetac_forget) * thetack_init(row, i);
					}
    					continue;
				}
				//加速度阶梯限幅
				//if (fabsf(acc_delta) > acc_change_limit) {
    				//	uk(0, i) = uk(1, i) + matrix::sign(acc_delta) * acc_change_limit;
				//}
				//保证加速度方向与速度方向一致
				//if(sign(uk(0,i))!=sign(_vel_sp(i))){
				//	uk(0,i)=0;
				//}
			}

   	     }
		else{//固定增益模式，调试用
		   	for(int i = 0;i < 3;i++){
				thetamk(0,i)=_mfac_vel_thetam(i);
		   	}
		   	for(int i = 0;i < 3;i++){
				for(int j = 0;j<3;j++){
				thetack(j,i)=_mfac_vel_thetac(j,i);}
		 	}
		  	 //求Hk
		  	 for (int j = 0; j < 3; j++) {
   		       	 	Hk(0,j) = -ek(0,j);
   	 	  	 }
   	 	  	 for (int j = 0; j < 3; j++) {
    		      	 	Hk(1,j)= ek(0,j) - ek(1,j);
    		  	 }
    		  	 for (int j = 0; j < 3; j++) {
    		      	 	Hk(2,j) = ek(1,j) - ek(2,j);
    		   	}
			for(int i=0;i<2;i++){
	 	  		matrix::Vector3f Hk_col = matrix::Vector3f(Hk.col(i));
	   			float uktemp = (get_thetack(thetack,i).transpose()*Hk_col)(0, 0);
	  	 		uk(0,i)=uk(1,i)+uktemp;
			}
		}
		//最后Z轴输出还是用PID
		Vector3f acc_sp_velocity = vel_error.emult(_gain_vel_p) + _vel_int - _vel_dot.emult(_gain_vel_d);
		//xy轴加速度限幅
		for(int i=0;i<2;i++){
			if(uk(0,i)>_mfac_vel_acc_xy_Limit||uk(0,i)<-_mfac_vel_acc_xy_Limit)
			uk(0,i)=sign(uk(0,i))*_mfac_vel_acc_xy_Limit;
		}
		//其余替换为MFAC计算输出量
		for(int i=0;i<2;i++){
			acc_sp_velocity(i)=uk(0,i);
		}
		ControlMath::addIfNotNanVector3f(_acc_sp, acc_sp_velocity);
	}
	//如果控制Z轴速度则多加一部分,否则用默认的PID控制

	print_counter++;
	if (print_counter >= 50) { // 250 Hz / 50 = 5 Hz
		print_counter=0;
    		PX4_INFO(
        		"thetacX=[%.2f %.2f %.2f], thetamX=%.3f",
        		(double)thetack(0,0),
        		(double)thetack(1,0),
        		(double)thetack(2,0),
        		(double)thetamk(0,0)
    		);
		PX4_INFO(
        		"thetacY=[%.2f %.2f %.2f], thetamY=%.3f",
        		(double)thetack(0,1),
        		(double)thetack(1,1),
        		(double)thetack(2,1),
        		(double)thetamk(0,1)
    		);
	}
	//Z轴速度限制
	_vel_int(2) = math::constrain(_vel_int(2), -CONSTANTS_ONE_G, CONSTANTS_ONE_G);
	//记录上一时刻速度
	for(int i=0;i<3;i++){
		_velk1(i)=_vel(i);
	}
	}
	//MFAC_EXIT:
	_accelerationControl();//计算推力

	// 垂直方向推力抗积分饱和
	if ((_thr_sp(2) >= -_lim_thr_min && vel_error(2) >= 0.f) ||
	    (_thr_sp(2) <= -_lim_thr_max && vel_error(2) <= 0.f)) {
		vel_error(2) = 0.f;
	}

	// Prioritize vertical control while keeping a horizontal margin
	const Vector2f thrust_sp_xy(_thr_sp);
	const float thrust_sp_xy_norm = thrust_sp_xy.norm();
	const float thrust_max_squared = math::sq(_lim_thr_max);

	// 垂直推力优先保障
	const float allocated_horizontal_thrust = math::min(thrust_sp_xy_norm, _lim_thr_xy_margin);
	const float thrust_z_max_squared = thrust_max_squared - math::sq(allocated_horizontal_thrust);

	// 垂直推力限幅
	_thr_sp(2) = math::max(_thr_sp(2), -sqrtf(thrust_z_max_squared));

	// 水平推力动态分配
	const float thrust_max_xy_squared = thrust_max_squared - math::sq(_thr_sp(2));
	float thrust_max_xy = 0.f;

	if (thrust_max_xy_squared > 0.f) {
		thrust_max_xy = sqrtf(thrust_max_xy_squared);
	}

	// Saturate thrust in horizontal direction
	if (thrust_sp_xy_norm > thrust_max_xy) {
		_thr_sp.xy() = thrust_sp_xy / thrust_sp_xy_norm * thrust_max_xy;
	}
	//抗积分饱和
	// Use tracking Anti-Windup for horizontal direction: during saturation, the integrator is used to unsaturate the output
	// see Anti-Reset Windup for PID controllers, L.Rundqwist, 1990
	const Vector2f acc_sp_xy_produced = Vector2f(_thr_sp) * (CONSTANTS_ONE_G / _hover_thrust);
	const float arw_gain = 2.f / _gain_vel_p(0);
	//当因为约束无法达到期望的速度时
	// The produced acceleration can be greater or smaller than the desired acceleration due to the saturations and the actual vertical thrust (computed independently).
	// The ARW loop needs to run if the signal is saturated only.
	const Vector2f acc_sp_xy = _acc_sp.xy();
	const Vector2f acc_limited_xy = (acc_sp_xy.norm_squared() > acc_sp_xy_produced.norm_squared())
					? acc_sp_xy_produced//饱和时使用实际能产生的加速度
					: acc_sp_xy;//没饱和就随意
	//使用期望加速度减去实际能产生的最大加速度乘以饱和增益来控制误差累积
	vel_error.xy() = Vector2f(vel_error) - arw_gain * (acc_sp_xy - acc_limited_xy);

	// Make sure integral doesn't get NAN
	ControlMath::setZeroIfNanVector3f(vel_error);
	// Update integral part of velocity control
	_vel_int += vel_error.emult(_gain_vel_i) * dt;
}






void PositionControl::_accelerationControl()//计算推力的
{
	// Assume standard acceleration due to gravity in vertical direction for attitude generation
	float z_specific_force = -CONSTANTS_ONE_G;

	if (!_decouple_horizontal_and_vertical_acceleration) {
		// Include vertical acceleration setpoint for better horizontal acceleration tracking
		z_specific_force += _acc_sp(2);
	}

	Vector3f body_z = Vector3f(-_acc_sp(0), -_acc_sp(1), -z_specific_force).normalized();
	ControlMath::limitTilt(body_z, Vector3f(0, 0, 1), _lim_tilt);
	// Convert to thrust assuming hover thrust produces standard gravity
	const float thrust_ned_z = _acc_sp(2) * (_hover_thrust / CONSTANTS_ONE_G) - _hover_thrust;


	// Project thrust to planned body attitude
	float cos_ned_body = (Vector3f(0, 0, 1).dot(body_z));

	// 防止出现 NaN / inf,给个最小值
    	if (!PX4_ISFINITE(cos_ned_body)) {
        	cos_ned_body = eps;
    	}

        float collective_thrust = math::min(thrust_ned_z / cos_ned_body, -_lim_thr_min);

	 if (cos_ned_body > eps) {
        // 正常情况：
        // 根据当前倾斜角，将期望的 NED 推力投影到机体方向
        	collective_thrust = thrust_ned_z / cos_ned_body;
    	} else {
        // 异常情况：
        // 当机体几乎水平或姿态异常时，避免除以接近 0 的值
        // 不做投影放大，使用保守的推力值
        	collective_thrust = thrust_ned_z;
    	}
	// 对总推力进行上下限约束（NED 坐标系下推力为负值）
    	// -_lim_thr_max ：最大向上推力
    	// -_lim_thr_min ：最小向上推力（接近 0）
    	const float thrust_min_allowed = -_lim_thr_max;
    	const float thrust_max_allowed = -_lim_thr_min;
    	collective_thrust = math::constrain(
       		collective_thrust,
        	thrust_min_allowed,
        	thrust_max_allowed
    	);

    	// 生成三轴推力设定值
	_thr_sp = body_z * collective_thrust;




	//print_counter++;
	//if (print_counter >= 250) { // 250 Hz / 50 = 5 Hz
    	//	print_counter = 0;
//
    	//	PX4_INFO(
        //		"accSp=[%.2f %.2f %.2f], thrust=%.3f",
        //		(double)_acc_sp(0),
        //		(double)_acc_sp(1),
        //		(double)_acc_sp(2),
        //		(double)_thr_sp(2)
    	//	);
	//}
}

bool PositionControl::_inputValid()
{
	bool valid = true;

	// Every axis x, y, z needs to have some setpoint
	for (int i = 0; i <= 2; i++) {
		valid = valid && (PX4_ISFINITE(_pos_sp(i)) || PX4_ISFINITE(_vel_sp(i)) || PX4_ISFINITE(_acc_sp(i)));
	}

	// x and y input setpoints always have to come in pairs
	valid = valid && (PX4_ISFINITE(_pos_sp(0)) == PX4_ISFINITE(_pos_sp(1)));
	valid = valid && (PX4_ISFINITE(_vel_sp(0)) == PX4_ISFINITE(_vel_sp(1)));
	valid = valid && (PX4_ISFINITE(_acc_sp(0)) == PX4_ISFINITE(_acc_sp(1)));

	// For each controlled state the estimate has to be valid
	for (int i = 0; i <= 2; i++) {
		if (PX4_ISFINITE(_pos_sp(i))) {
			valid = valid && PX4_ISFINITE(_pos(i));
		}

		if (PX4_ISFINITE(_vel_sp(i))) {
			valid = valid && PX4_ISFINITE(_vel(i)) && PX4_ISFINITE(_vel_dot(i));
		}
	}

	return valid;
}

void PositionControl::getLocalPositionSetpoint(vehicle_local_position_setpoint_s &local_position_setpoint) const
{
	local_position_setpoint.x = _pos_sp(0);
	local_position_setpoint.y = _pos_sp(1);
	local_position_setpoint.z = _pos_sp(2);
	local_position_setpoint.yaw = _yaw_sp;
	local_position_setpoint.yawspeed = _yawspeed_sp;
	local_position_setpoint.vx = _vel_sp(0);
	local_position_setpoint.vy = _vel_sp(1);
	local_position_setpoint.vz = _vel_sp(2);
	_acc_sp.copyTo(local_position_setpoint.acceleration);
	_thr_sp.copyTo(local_position_setpoint.thrust);
}

void PositionControl::getAttitudeSetpoint(vehicle_attitude_setpoint_s &attitude_setpoint) const
{
	ControlMath::thrustToAttitude(_thr_sp, _yaw_sp, attitude_setpoint);
	attitude_setpoint.yaw_sp_move_rate = _yawspeed_sp;
}
