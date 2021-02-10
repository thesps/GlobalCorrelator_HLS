-- ==============================================================
-- RTL generated by Vivado(TM) HLS - High-Level Synthesis from C, C++ and OpenCL
-- Version: 2019.2
-- Copyright (C) 1986-2019 Xilinx, Inc. All Rights Reserved.
-- 
-- ===========================================================

library IEEE;
use IEEE.std_logic_1164.all;
use IEEE.numeric_std.all;

entity reduce_14 is
port (
    ap_clk : IN STD_LOGIC;
    ap_rst : IN STD_LOGIC;
    x_0_V_read : IN STD_LOGIC_VECTOR (21 downto 0);
    x_1_V_read : IN STD_LOGIC_VECTOR (21 downto 0);
    x_2_V_read : IN STD_LOGIC_VECTOR (21 downto 0);
    x_3_V_read : IN STD_LOGIC_VECTOR (21 downto 0);
    x_4_V_read : IN STD_LOGIC_VECTOR (21 downto 0);
    x_5_V_read : IN STD_LOGIC_VECTOR (21 downto 0);
    x_6_V_read : IN STD_LOGIC_VECTOR (21 downto 0);
    x_7_V_read : IN STD_LOGIC_VECTOR (21 downto 0);
    ap_return : OUT STD_LOGIC_VECTOR (21 downto 0);
    ap_ce : IN STD_LOGIC );
end;


architecture behav of reduce_14 is 
    constant ap_const_logic_1 : STD_LOGIC := '1';
    constant ap_const_boolean_1 : BOOLEAN := true;
    constant ap_const_boolean_0 : BOOLEAN := false;
    constant ap_const_logic_0 : STD_LOGIC := '0';

    signal p_Val2_s_reduce_16_fu_76_ap_return : STD_LOGIC_VECTOR (21 downto 0);
    signal p_Val2_s_reg_104 : STD_LOGIC_VECTOR (21 downto 0);
    signal ap_block_state1_pp0_stage0_iter0 : BOOLEAN;
    signal ap_block_state2_pp0_stage0_iter1 : BOOLEAN;
    signal ap_block_pp0_stage0_11001 : BOOLEAN;
    signal p_Val2_7_reduce_16_fu_88_ap_return : STD_LOGIC_VECTOR (21 downto 0);
    signal p_Val2_7_reg_109 : STD_LOGIC_VECTOR (21 downto 0);
    signal p_Val2_s_reduce_16_fu_76_ap_ready : STD_LOGIC;
    signal p_Val2_7_reduce_16_fu_88_ap_ready : STD_LOGIC;
    signal ap_block_pp0_stage0 : BOOLEAN;

    component reduce_16 IS
    port (
        ap_ready : OUT STD_LOGIC;
        x_0_V_read : IN STD_LOGIC_VECTOR (21 downto 0);
        x_1_V_read : IN STD_LOGIC_VECTOR (21 downto 0);
        x_2_V_read : IN STD_LOGIC_VECTOR (21 downto 0);
        x_3_V_read : IN STD_LOGIC_VECTOR (21 downto 0);
        ap_return : OUT STD_LOGIC_VECTOR (21 downto 0) );
    end component;



begin
    p_Val2_s_reduce_16_fu_76 : component reduce_16
    port map (
        ap_ready => p_Val2_s_reduce_16_fu_76_ap_ready,
        x_0_V_read => x_0_V_read,
        x_1_V_read => x_1_V_read,
        x_2_V_read => x_2_V_read,
        x_3_V_read => x_3_V_read,
        ap_return => p_Val2_s_reduce_16_fu_76_ap_return);

    p_Val2_7_reduce_16_fu_88 : component reduce_16
    port map (
        ap_ready => p_Val2_7_reduce_16_fu_88_ap_ready,
        x_0_V_read => x_4_V_read,
        x_1_V_read => x_5_V_read,
        x_2_V_read => x_6_V_read,
        x_3_V_read => x_7_V_read,
        ap_return => p_Val2_7_reduce_16_fu_88_ap_return);




    process (ap_clk)
    begin
        if (ap_clk'event and ap_clk = '1') then
            if (((ap_const_boolean_0 = ap_block_pp0_stage0_11001) and (ap_const_logic_1 = ap_ce))) then
                p_Val2_7_reg_109 <= p_Val2_7_reduce_16_fu_88_ap_return;
                p_Val2_s_reg_104 <= p_Val2_s_reduce_16_fu_76_ap_return;
            end if;
        end if;
    end process;
        ap_block_pp0_stage0 <= not((ap_const_boolean_1 = ap_const_boolean_1));
        ap_block_pp0_stage0_11001 <= not((ap_const_boolean_1 = ap_const_boolean_1));
        ap_block_state1_pp0_stage0_iter0 <= not((ap_const_boolean_1 = ap_const_boolean_1));
        ap_block_state2_pp0_stage0_iter1 <= not((ap_const_boolean_1 = ap_const_boolean_1));
    ap_return <= std_logic_vector(unsigned(p_Val2_7_reg_109) + unsigned(p_Val2_s_reg_104));
end behav;
