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
    ap_clk : IN STD_LOGIC;
    ap_rst : IN STD_LOGIC;
    x_0_hwPt_V_read : IN STD_LOGIC_VECTOR (15 downto 0);
    x_1_hwPt_V_read : IN STD_LOGIC_VECTOR (15 downto 0);
    x_2_hwPt_V_read : IN STD_LOGIC_VECTOR (15 downto 0);
    x_3_hwPt_V_read : IN STD_LOGIC_VECTOR (15 downto 0);
    x_0_hwEta_V_read : IN STD_LOGIC_VECTOR (9 downto 0);
    x_1_hwEta_V_read : IN STD_LOGIC_VECTOR (9 downto 0);
    x_2_hwEta_V_read : IN STD_LOGIC_VECTOR (9 downto 0);
    x_3_hwEta_V_read : IN STD_LOGIC_VECTOR (9 downto 0);
    x_0_hwPhi_V_read : IN STD_LOGIC_VECTOR (9 downto 0);
    x_1_hwPhi_V_read : IN STD_LOGIC_VECTOR (9 downto 0);
    x_2_hwPhi_V_read : IN STD_LOGIC_VECTOR (9 downto 0);
    x_3_hwPhi_V_read : IN STD_LOGIC_VECTOR (9 downto 0);
    ap_return_0 : OUT STD_LOGIC_VECTOR (15 downto 0);
    ap_return_1 : OUT STD_LOGIC_VECTOR (9 downto 0);
    ap_return_2 : OUT STD_LOGIC_VECTOR (9 downto 0);
    ap_ce : IN STD_LOGIC );
end;


architecture behav of reduce_2 is 
    constant ap_const_logic_1 : STD_LOGIC := '1';
    constant ap_const_logic_0 : STD_LOGIC := '0';
    constant ap_const_boolean_1 : BOOLEAN := true;
    constant ap_const_boolean_0 : BOOLEAN := false;
    constant ap_const_lv1_1 : STD_LOGIC_VECTOR (0 downto 0) := "1";

    signal call_ret5_reduce_4_fu_114_ap_return_0 : STD_LOGIC_VECTOR (15 downto 0);
    signal call_ret5_reduce_4_fu_114_ap_return_1 : STD_LOGIC_VECTOR (9 downto 0);
    signal call_ret5_reduce_4_fu_114_ap_return_2 : STD_LOGIC_VECTOR (9 downto 0);
    signal call_ret5_reg_217_1 : STD_LOGIC_VECTOR (9 downto 0);
    signal call_ret5_reg_217_2 : STD_LOGIC_VECTOR (9 downto 0);
    signal ap_block_state1_pp0_stage0_iter0 : BOOLEAN;
    signal ap_block_state2_pp0_stage0_iter1 : BOOLEAN;
    signal ap_block_pp0_stage0_11001 : BOOLEAN;
    signal call_ret6_reduce_4_fu_130_ap_return_0 : STD_LOGIC_VECTOR (15 downto 0);
    signal call_ret6_reduce_4_fu_130_ap_return_1 : STD_LOGIC_VECTOR (9 downto 0);
    signal call_ret6_reduce_4_fu_130_ap_return_2 : STD_LOGIC_VECTOR (9 downto 0);
    signal call_ret6_reg_223_1 : STD_LOGIC_VECTOR (9 downto 0);
    signal call_ret6_reg_223_2 : STD_LOGIC_VECTOR (9 downto 0);
    signal xor_ln1496_fu_160_p2 : STD_LOGIC_VECTOR (0 downto 0);
    signal xor_ln1496_reg_229 : STD_LOGIC_VECTOR (0 downto 0);
    signal select_ln93_fu_166_p3 : STD_LOGIC_VECTOR (15 downto 0);
    signal select_ln93_reg_235 : STD_LOGIC_VECTOR (15 downto 0);
    signal call_ret5_reduce_4_fu_114_ap_ready : STD_LOGIC;
    signal call_ret6_reduce_4_fu_130_ap_ready : STD_LOGIC;
    signal ap_block_pp0_stage0 : BOOLEAN;
    signal icmp_ln1496_fu_154_p2 : STD_LOGIC_VECTOR (0 downto 0);
    signal select_ln93_5_fu_186_p3 : STD_LOGIC_VECTOR (9 downto 0);
    signal select_ln93_6_fu_193_p3 : STD_LOGIC_VECTOR (9 downto 0);
    signal ap_ce_reg : STD_LOGIC;
    signal ap_return_0_int_reg : STD_LOGIC_VECTOR (15 downto 0);
    signal ap_return_1_int_reg : STD_LOGIC_VECTOR (9 downto 0);
    signal ap_return_2_int_reg : STD_LOGIC_VECTOR (9 downto 0);

    component reduce_4 IS
    port (
        ap_ready : OUT STD_LOGIC;
        x_0_hwPt_V_read : IN STD_LOGIC_VECTOR (15 downto 0);
        x_1_hwPt_V_read : IN STD_LOGIC_VECTOR (15 downto 0);
        x_0_hwEta_V_read : IN STD_LOGIC_VECTOR (9 downto 0);
        x_1_hwEta_V_read : IN STD_LOGIC_VECTOR (9 downto 0);
        x_0_hwPhi_V_read : IN STD_LOGIC_VECTOR (9 downto 0);
        x_1_hwPhi_V_read : IN STD_LOGIC_VECTOR (9 downto 0);
        ap_return_0 : OUT STD_LOGIC_VECTOR (15 downto 0);
        ap_return_1 : OUT STD_LOGIC_VECTOR (9 downto 0);
        ap_return_2 : OUT STD_LOGIC_VECTOR (9 downto 0) );
    end component;



begin
    call_ret5_reduce_4_fu_114 : component reduce_4
    port map (
        ap_ready => call_ret5_reduce_4_fu_114_ap_ready,
        x_0_hwPt_V_read => x_0_hwPt_V_read,
        x_1_hwPt_V_read => x_1_hwPt_V_read,
        x_0_hwEta_V_read => x_0_hwEta_V_read,
        x_1_hwEta_V_read => x_1_hwEta_V_read,
        x_0_hwPhi_V_read => x_0_hwPhi_V_read,
        x_1_hwPhi_V_read => x_1_hwPhi_V_read,
        ap_return_0 => call_ret5_reduce_4_fu_114_ap_return_0,
        ap_return_1 => call_ret5_reduce_4_fu_114_ap_return_1,
        ap_return_2 => call_ret5_reduce_4_fu_114_ap_return_2);

    call_ret6_reduce_4_fu_130 : component reduce_4
    port map (
        ap_ready => call_ret6_reduce_4_fu_130_ap_ready,
        x_0_hwPt_V_read => x_2_hwPt_V_read,
        x_1_hwPt_V_read => x_3_hwPt_V_read,
        x_0_hwEta_V_read => x_2_hwEta_V_read,
        x_1_hwEta_V_read => x_3_hwEta_V_read,
        x_0_hwPhi_V_read => x_2_hwPhi_V_read,
        x_1_hwPhi_V_read => x_3_hwPhi_V_read,
        ap_return_0 => call_ret6_reduce_4_fu_130_ap_return_0,
        ap_return_1 => call_ret6_reduce_4_fu_130_ap_return_1,
        ap_return_2 => call_ret6_reduce_4_fu_130_ap_return_2);





    ap_ce_reg_assign_proc : process (ap_clk)
    begin
        if (ap_clk'event and ap_clk = '1') then
            ap_ce_reg <= ap_ce;
        end if;
    end process;
    process (ap_clk)
    begin
        if (ap_clk'event and ap_clk = '1') then
            if ((ap_const_logic_1 = ap_ce_reg)) then
                ap_return_0_int_reg <= select_ln93_reg_235;
                ap_return_1_int_reg <= select_ln93_5_fu_186_p3;
                ap_return_2_int_reg <= select_ln93_6_fu_193_p3;
            end if;
        end if;
    end process;
    process (ap_clk)
    begin
        if (ap_clk'event and ap_clk = '1') then
            if (((ap_const_boolean_0 = ap_block_pp0_stage0_11001) and (ap_const_logic_1 = ap_const_logic_1))) then
                call_ret5_reg_217_1 <= call_ret5_reduce_4_fu_114_ap_return_1;
                call_ret5_reg_217_2 <= call_ret5_reduce_4_fu_114_ap_return_2;
                call_ret6_reg_223_1 <= call_ret6_reduce_4_fu_130_ap_return_1;
                call_ret6_reg_223_2 <= call_ret6_reduce_4_fu_130_ap_return_2;
                select_ln93_reg_235 <= select_ln93_fu_166_p3;
                xor_ln1496_reg_229 <= xor_ln1496_fu_160_p2;
            end if;
        end if;
    end process;
        ap_block_pp0_stage0 <= not((ap_const_boolean_1 = ap_const_boolean_1));
        ap_block_pp0_stage0_11001 <= not((ap_const_boolean_1 = ap_const_boolean_1));
        ap_block_state1_pp0_stage0_iter0 <= not((ap_const_boolean_1 = ap_const_boolean_1));
        ap_block_state2_pp0_stage0_iter1 <= not((ap_const_boolean_1 = ap_const_boolean_1));

    ap_return_0_assign_proc : process(select_ln93_reg_235, ap_ce_reg, ap_return_0_int_reg)
    begin
        if ((ap_const_logic_0 = ap_ce_reg)) then 
            ap_return_0 <= ap_return_0_int_reg;
        elsif ((ap_const_logic_1 = ap_ce_reg)) then 
            ap_return_0 <= select_ln93_reg_235;
        end if; 
    end process;


    ap_return_1_assign_proc : process(select_ln93_5_fu_186_p3, ap_ce_reg, ap_return_1_int_reg)
    begin
        if ((ap_const_logic_0 = ap_ce_reg)) then 
            ap_return_1 <= ap_return_1_int_reg;
        elsif ((ap_const_logic_1 = ap_ce_reg)) then 
            ap_return_1 <= select_ln93_5_fu_186_p3;
        end if; 
    end process;


    ap_return_2_assign_proc : process(select_ln93_6_fu_193_p3, ap_ce_reg, ap_return_2_int_reg)
    begin
        if ((ap_const_logic_0 = ap_ce_reg)) then 
            ap_return_2 <= ap_return_2_int_reg;
        elsif ((ap_const_logic_1 = ap_ce_reg)) then 
            ap_return_2 <= select_ln93_6_fu_193_p3;
        end if; 
    end process;

    icmp_ln1496_fu_154_p2 <= "1" when (unsigned(call_ret5_reduce_4_fu_114_ap_return_0) < unsigned(call_ret6_reduce_4_fu_130_ap_return_0)) else "0";
    select_ln93_5_fu_186_p3 <= 
        call_ret5_reg_217_1 when (xor_ln1496_reg_229(0) = '1') else 
        call_ret6_reg_223_1;
    select_ln93_6_fu_193_p3 <= 
        call_ret5_reg_217_2 when (xor_ln1496_reg_229(0) = '1') else 
        call_ret6_reg_223_2;
    select_ln93_fu_166_p3 <= 
        call_ret5_reduce_4_fu_114_ap_return_0 when (xor_ln1496_fu_160_p2(0) = '1') else 
        call_ret6_reduce_4_fu_130_ap_return_0;
    xor_ln1496_fu_160_p2 <= (icmp_ln1496_fu_154_p2 xor ap_const_lv1_1);
end behav;
