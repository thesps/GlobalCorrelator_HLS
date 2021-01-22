-- ==============================================================
-- RTL generated by Vivado(TM) HLS - High-Level Synthesis from C, C++ and OpenCL
-- Version: 2019.2
-- Copyright (C) 1986-2019 Xilinx, Inc. All Rights Reserved.
-- 
-- ===========================================================

library IEEE;
use IEEE.std_logic_1164.all;
use IEEE.numeric_std.all;

entity reduce_2 is
port (
    ap_ready : OUT STD_LOGIC;
    x_0_V_read : IN STD_LOGIC_VECTOR (0 downto 0);
    x_1_V_read : IN STD_LOGIC_VECTOR (0 downto 0);
    x_2_V_read : IN STD_LOGIC_VECTOR (0 downto 0);
    x_3_V_read : IN STD_LOGIC_VECTOR (0 downto 0);
    ap_return : OUT STD_LOGIC_VECTOR (2 downto 0) );
end;


architecture behav of reduce_2 is 
    constant ap_const_logic_1 : STD_LOGIC := '1';
    constant ap_const_boolean_1 : BOOLEAN := true;
    constant ap_const_logic_0 : STD_LOGIC := '0';

    signal p_s_reduce_4_fu_44_ap_ready : STD_LOGIC;
    signal p_s_reduce_4_fu_44_ap_return : STD_LOGIC_VECTOR (1 downto 0);
    signal p_3_reduce_4_fu_52_ap_ready : STD_LOGIC;
    signal p_3_reduce_4_fu_52_ap_return : STD_LOGIC_VECTOR (1 downto 0);
    signal zext_ln209_fu_60_p1 : STD_LOGIC_VECTOR (2 downto 0);
    signal zext_ln209_1_fu_64_p1 : STD_LOGIC_VECTOR (2 downto 0);

    component reduce_4 IS
    port (
        ap_ready : OUT STD_LOGIC;
        x_0_V_read : IN STD_LOGIC_VECTOR (0 downto 0);
        x_1_V_read : IN STD_LOGIC_VECTOR (0 downto 0);
        ap_return : OUT STD_LOGIC_VECTOR (1 downto 0) );
    end component;



begin
    p_s_reduce_4_fu_44 : component reduce_4
    port map (
        ap_ready => p_s_reduce_4_fu_44_ap_ready,
        x_0_V_read => x_0_V_read,
        x_1_V_read => x_1_V_read,
        ap_return => p_s_reduce_4_fu_44_ap_return);

    p_3_reduce_4_fu_52 : component reduce_4
    port map (
        ap_ready => p_3_reduce_4_fu_52_ap_ready,
        x_0_V_read => x_2_V_read,
        x_1_V_read => x_3_V_read,
        ap_return => p_3_reduce_4_fu_52_ap_return);




    ap_ready <= ap_const_logic_1;
    ap_return <= std_logic_vector(unsigned(zext_ln209_fu_60_p1) + unsigned(zext_ln209_1_fu_64_p1));
    zext_ln209_1_fu_64_p1 <= std_logic_vector(IEEE.numeric_std.resize(unsigned(p_3_reduce_4_fu_52_ap_return),3));
    zext_ln209_fu_60_p1 <= std_logic_vector(IEEE.numeric_std.resize(unsigned(p_s_reduce_4_fu_44_ap_return),3));
end behav;
