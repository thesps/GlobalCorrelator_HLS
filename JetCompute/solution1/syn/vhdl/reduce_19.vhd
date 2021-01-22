-- ==============================================================
-- RTL generated by Vivado(TM) HLS - High-Level Synthesis from C, C++ and OpenCL
-- Version: 2019.2
-- Copyright (C) 1986-2019 Xilinx, Inc. All Rights Reserved.
-- 
-- ===========================================================

library IEEE;
use IEEE.std_logic_1164.all;
use IEEE.numeric_std.all;

entity reduce_19 is
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
    x_8_V_read : IN STD_LOGIC_VECTOR (21 downto 0);
    x_9_V_read : IN STD_LOGIC_VECTOR (21 downto 0);
    x_10_V_read : IN STD_LOGIC_VECTOR (21 downto 0);
    x_11_V_read : IN STD_LOGIC_VECTOR (21 downto 0);
    x_12_V_read : IN STD_LOGIC_VECTOR (21 downto 0);
    x_13_V_read : IN STD_LOGIC_VECTOR (21 downto 0);
    x_14_V_read : IN STD_LOGIC_VECTOR (21 downto 0);
    x_15_V_read : IN STD_LOGIC_VECTOR (21 downto 0);
    ap_return : OUT STD_LOGIC_VECTOR (21 downto 0);
    ap_ce : IN STD_LOGIC );
end;


architecture behav of reduce_19 is 
    constant ap_const_logic_1 : STD_LOGIC := '1';
    constant ap_const_logic_0 : STD_LOGIC := '0';
    constant ap_const_boolean_1 : BOOLEAN := true;
    constant ap_const_boolean_0 : BOOLEAN := false;

    signal p_Val2_s_reduce_14_fu_140_ap_return : STD_LOGIC_VECTOR (21 downto 0);
    signal p_Val2_s_reg_184 : STD_LOGIC_VECTOR (21 downto 0);
    signal ap_block_state1_pp0_stage0_iter0 : BOOLEAN;
    signal ap_block_state2_pp0_stage0_iter1 : BOOLEAN;
    signal ap_block_pp0_stage0_11001 : BOOLEAN;
    signal p_Val2_4_reduce_14_fu_160_ap_return : STD_LOGIC_VECTOR (21 downto 0);
    signal p_Val2_4_reg_189 : STD_LOGIC_VECTOR (21 downto 0);
    signal p_Val2_s_reduce_14_fu_140_ap_ready : STD_LOGIC;
    signal p_Val2_4_reduce_14_fu_160_ap_ready : STD_LOGIC;
    signal ap_block_pp0_stage0 : BOOLEAN;
    signal add_ln703_fu_180_p2 : STD_LOGIC_VECTOR (21 downto 0);
    signal ap_ce_reg : STD_LOGIC;
    signal ap_return_int_reg : STD_LOGIC_VECTOR (21 downto 0);

    component reduce_14 IS
    port (
        ap_ready : OUT STD_LOGIC;
        x_0_V_read : IN STD_LOGIC_VECTOR (21 downto 0);
        x_1_V_read : IN STD_LOGIC_VECTOR (21 downto 0);
        x_2_V_read : IN STD_LOGIC_VECTOR (21 downto 0);
        x_3_V_read : IN STD_LOGIC_VECTOR (21 downto 0);
        x_4_V_read : IN STD_LOGIC_VECTOR (21 downto 0);
        x_5_V_read : IN STD_LOGIC_VECTOR (21 downto 0);
        x_6_V_read : IN STD_LOGIC_VECTOR (21 downto 0);
        x_7_V_read : IN STD_LOGIC_VECTOR (21 downto 0);
        ap_return : OUT STD_LOGIC_VECTOR (21 downto 0) );
    end component;



begin
    p_Val2_s_reduce_14_fu_140 : component reduce_14
    port map (
        ap_ready => p_Val2_s_reduce_14_fu_140_ap_ready,
        x_0_V_read => x_0_V_read,
        x_1_V_read => x_1_V_read,
        x_2_V_read => x_2_V_read,
        x_3_V_read => x_3_V_read,
        x_4_V_read => x_4_V_read,
        x_5_V_read => x_5_V_read,
        x_6_V_read => x_6_V_read,
        x_7_V_read => x_7_V_read,
        ap_return => p_Val2_s_reduce_14_fu_140_ap_return);

    p_Val2_4_reduce_14_fu_160 : component reduce_14
    port map (
        ap_ready => p_Val2_4_reduce_14_fu_160_ap_ready,
        x_0_V_read => x_8_V_read,
        x_1_V_read => x_9_V_read,
        x_2_V_read => x_10_V_read,
        x_3_V_read => x_11_V_read,
        x_4_V_read => x_12_V_read,
        x_5_V_read => x_13_V_read,
        x_6_V_read => x_14_V_read,
        x_7_V_read => x_15_V_read,
        ap_return => p_Val2_4_reduce_14_fu_160_ap_return);





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
                ap_return_int_reg <= add_ln703_fu_180_p2;
            end if;
        end if;
    end process;
    process (ap_clk)
    begin
        if (ap_clk'event and ap_clk = '1') then
            if (((ap_const_logic_1 = ap_const_logic_1) and (ap_const_boolean_0 = ap_block_pp0_stage0_11001))) then
                p_Val2_4_reg_189 <= p_Val2_4_reduce_14_fu_160_ap_return;
                p_Val2_s_reg_184 <= p_Val2_s_reduce_14_fu_140_ap_return;
            end if;
        end if;
    end process;
    add_ln703_fu_180_p2 <= std_logic_vector(unsigned(p_Val2_4_reg_189) + unsigned(p_Val2_s_reg_184));
        ap_block_pp0_stage0 <= not((ap_const_boolean_1 = ap_const_boolean_1));
        ap_block_pp0_stage0_11001 <= not((ap_const_boolean_1 = ap_const_boolean_1));
        ap_block_state1_pp0_stage0_iter0 <= not((ap_const_boolean_1 = ap_const_boolean_1));
        ap_block_state2_pp0_stage0_iter1 <= not((ap_const_boolean_1 = ap_const_boolean_1));

    ap_return_assign_proc : process(add_ln703_fu_180_p2, ap_ce_reg, ap_return_int_reg)
    begin
        if ((ap_const_logic_0 = ap_ce_reg)) then 
            ap_return <= ap_return_int_reg;
        elsif ((ap_const_logic_1 = ap_ce_reg)) then 
            ap_return <= add_ln703_fu_180_p2;
        end if; 
    end process;

end behav;
