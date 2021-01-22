-- ==============================================================
-- RTL generated by Vivado(TM) HLS - High-Level Synthesis from C, C++ and OpenCL
-- Version: 2019.2
-- Copyright (C) 1986-2019 Xilinx, Inc. All Rights Reserved.
-- 
-- ===========================================================

library IEEE;
use IEEE.std_logic_1164.all;
use IEEE.numeric_std.all;

entity reduce_15 is
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
    x_16_V_read : IN STD_LOGIC_VECTOR (21 downto 0);
    x_17_V_read : IN STD_LOGIC_VECTOR (21 downto 0);
    x_18_V_read : IN STD_LOGIC_VECTOR (21 downto 0);
    x_19_V_read : IN STD_LOGIC_VECTOR (21 downto 0);
    x_20_V_read : IN STD_LOGIC_VECTOR (21 downto 0);
    x_21_V_read : IN STD_LOGIC_VECTOR (21 downto 0);
    x_22_V_read : IN STD_LOGIC_VECTOR (21 downto 0);
    x_23_V_read : IN STD_LOGIC_VECTOR (21 downto 0);
    x_24_V_read : IN STD_LOGIC_VECTOR (21 downto 0);
    x_25_V_read : IN STD_LOGIC_VECTOR (21 downto 0);
    x_26_V_read : IN STD_LOGIC_VECTOR (21 downto 0);
    x_27_V_read : IN STD_LOGIC_VECTOR (21 downto 0);
    x_28_V_read : IN STD_LOGIC_VECTOR (21 downto 0);
    x_29_V_read : IN STD_LOGIC_VECTOR (21 downto 0);
    x_30_V_read : IN STD_LOGIC_VECTOR (21 downto 0);
    x_31_V_read : IN STD_LOGIC_VECTOR (21 downto 0);
    x_32_V_read : IN STD_LOGIC_VECTOR (21 downto 0);
    x_33_V_read : IN STD_LOGIC_VECTOR (21 downto 0);
    x_34_V_read : IN STD_LOGIC_VECTOR (21 downto 0);
    x_35_V_read : IN STD_LOGIC_VECTOR (21 downto 0);
    x_36_V_read : IN STD_LOGIC_VECTOR (21 downto 0);
    x_37_V_read : IN STD_LOGIC_VECTOR (21 downto 0);
    x_38_V_read : IN STD_LOGIC_VECTOR (21 downto 0);
    x_39_V_read : IN STD_LOGIC_VECTOR (21 downto 0);
    x_40_V_read : IN STD_LOGIC_VECTOR (21 downto 0);
    x_41_V_read : IN STD_LOGIC_VECTOR (21 downto 0);
    x_42_V_read : IN STD_LOGIC_VECTOR (21 downto 0);
    x_43_V_read : IN STD_LOGIC_VECTOR (21 downto 0);
    x_44_V_read : IN STD_LOGIC_VECTOR (21 downto 0);
    x_45_V_read : IN STD_LOGIC_VECTOR (21 downto 0);
    x_46_V_read : IN STD_LOGIC_VECTOR (21 downto 0);
    x_47_V_read : IN STD_LOGIC_VECTOR (21 downto 0);
    x_48_V_read : IN STD_LOGIC_VECTOR (21 downto 0);
    x_49_V_read : IN STD_LOGIC_VECTOR (21 downto 0);
    x_50_V_read : IN STD_LOGIC_VECTOR (21 downto 0);
    x_51_V_read : IN STD_LOGIC_VECTOR (21 downto 0);
    x_52_V_read : IN STD_LOGIC_VECTOR (21 downto 0);
    x_53_V_read : IN STD_LOGIC_VECTOR (21 downto 0);
    x_54_V_read : IN STD_LOGIC_VECTOR (21 downto 0);
    x_55_V_read : IN STD_LOGIC_VECTOR (21 downto 0);
    x_56_V_read : IN STD_LOGIC_VECTOR (21 downto 0);
    x_57_V_read : IN STD_LOGIC_VECTOR (21 downto 0);
    x_58_V_read : IN STD_LOGIC_VECTOR (21 downto 0);
    x_59_V_read : IN STD_LOGIC_VECTOR (21 downto 0);
    x_60_V_read : IN STD_LOGIC_VECTOR (21 downto 0);
    x_61_V_read : IN STD_LOGIC_VECTOR (21 downto 0);
    x_62_V_read : IN STD_LOGIC_VECTOR (21 downto 0);
    x_63_V_read : IN STD_LOGIC_VECTOR (21 downto 0);
    ap_return : OUT STD_LOGIC_VECTOR (21 downto 0);
    ap_ce : IN STD_LOGIC );
end;


architecture behav of reduce_15 is 
    constant ap_const_logic_1 : STD_LOGIC := '1';
    constant ap_const_logic_0 : STD_LOGIC := '0';
    constant ap_const_boolean_1 : BOOLEAN := true;
    constant ap_const_boolean_0 : BOOLEAN := false;

    signal grp_reduce_17_fu_524_ap_return : STD_LOGIC_VECTOR (21 downto 0);
    signal grp_reduce_17_fu_524_ap_ce : STD_LOGIC;
    signal ap_block_state1_pp0_stage0_iter0_ignore_call65 : BOOLEAN;
    signal ap_block_state2_pp0_stage0_iter1_ignore_call65 : BOOLEAN;
    signal ap_block_pp0_stage0_11001_ignoreCallOp67 : BOOLEAN;
    signal grp_reduce_17_fu_592_ap_return : STD_LOGIC_VECTOR (21 downto 0);
    signal grp_reduce_17_fu_592_ap_ce : STD_LOGIC;
    signal ap_block_state1_pp0_stage0_iter0_ignore_call66 : BOOLEAN;
    signal ap_block_state2_pp0_stage0_iter1_ignore_call66 : BOOLEAN;
    signal ap_block_pp0_stage0_11001_ignoreCallOp68 : BOOLEAN;
    signal ap_block_pp0_stage0 : BOOLEAN;
    signal ap_block_state1_pp0_stage0_iter0 : BOOLEAN;
    signal ap_block_state2_pp0_stage0_iter1 : BOOLEAN;
    signal ap_block_pp0_stage0_11001 : BOOLEAN;
    signal add_ln703_fu_660_p2 : STD_LOGIC_VECTOR (21 downto 0);
    signal ap_ce_reg : STD_LOGIC;
    signal ap_return_int_reg : STD_LOGIC_VECTOR (21 downto 0);

    component reduce_17 IS
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
        x_16_V_read : IN STD_LOGIC_VECTOR (21 downto 0);
        x_17_V_read : IN STD_LOGIC_VECTOR (21 downto 0);
        x_18_V_read : IN STD_LOGIC_VECTOR (21 downto 0);
        x_19_V_read : IN STD_LOGIC_VECTOR (21 downto 0);
        x_20_V_read : IN STD_LOGIC_VECTOR (21 downto 0);
        x_21_V_read : IN STD_LOGIC_VECTOR (21 downto 0);
        x_22_V_read : IN STD_LOGIC_VECTOR (21 downto 0);
        x_23_V_read : IN STD_LOGIC_VECTOR (21 downto 0);
        x_24_V_read : IN STD_LOGIC_VECTOR (21 downto 0);
        x_25_V_read : IN STD_LOGIC_VECTOR (21 downto 0);
        x_26_V_read : IN STD_LOGIC_VECTOR (21 downto 0);
        x_27_V_read : IN STD_LOGIC_VECTOR (21 downto 0);
        x_28_V_read : IN STD_LOGIC_VECTOR (21 downto 0);
        x_29_V_read : IN STD_LOGIC_VECTOR (21 downto 0);
        x_30_V_read : IN STD_LOGIC_VECTOR (21 downto 0);
        x_31_V_read : IN STD_LOGIC_VECTOR (21 downto 0);
        ap_return : OUT STD_LOGIC_VECTOR (21 downto 0);
        ap_ce : IN STD_LOGIC );
    end component;



begin
    grp_reduce_17_fu_524 : component reduce_17
    port map (
        ap_clk => ap_clk,
        ap_rst => ap_rst,
        x_0_V_read => x_0_V_read,
        x_1_V_read => x_1_V_read,
        x_2_V_read => x_2_V_read,
        x_3_V_read => x_3_V_read,
        x_4_V_read => x_4_V_read,
        x_5_V_read => x_5_V_read,
        x_6_V_read => x_6_V_read,
        x_7_V_read => x_7_V_read,
        x_8_V_read => x_8_V_read,
        x_9_V_read => x_9_V_read,
        x_10_V_read => x_10_V_read,
        x_11_V_read => x_11_V_read,
        x_12_V_read => x_12_V_read,
        x_13_V_read => x_13_V_read,
        x_14_V_read => x_14_V_read,
        x_15_V_read => x_15_V_read,
        x_16_V_read => x_16_V_read,
        x_17_V_read => x_17_V_read,
        x_18_V_read => x_18_V_read,
        x_19_V_read => x_19_V_read,
        x_20_V_read => x_20_V_read,
        x_21_V_read => x_21_V_read,
        x_22_V_read => x_22_V_read,
        x_23_V_read => x_23_V_read,
        x_24_V_read => x_24_V_read,
        x_25_V_read => x_25_V_read,
        x_26_V_read => x_26_V_read,
        x_27_V_read => x_27_V_read,
        x_28_V_read => x_28_V_read,
        x_29_V_read => x_29_V_read,
        x_30_V_read => x_30_V_read,
        x_31_V_read => x_31_V_read,
        ap_return => grp_reduce_17_fu_524_ap_return,
        ap_ce => grp_reduce_17_fu_524_ap_ce);

    grp_reduce_17_fu_592 : component reduce_17
    port map (
        ap_clk => ap_clk,
        ap_rst => ap_rst,
        x_0_V_read => x_32_V_read,
        x_1_V_read => x_33_V_read,
        x_2_V_read => x_34_V_read,
        x_3_V_read => x_35_V_read,
        x_4_V_read => x_36_V_read,
        x_5_V_read => x_37_V_read,
        x_6_V_read => x_38_V_read,
        x_7_V_read => x_39_V_read,
        x_8_V_read => x_40_V_read,
        x_9_V_read => x_41_V_read,
        x_10_V_read => x_42_V_read,
        x_11_V_read => x_43_V_read,
        x_12_V_read => x_44_V_read,
        x_13_V_read => x_45_V_read,
        x_14_V_read => x_46_V_read,
        x_15_V_read => x_47_V_read,
        x_16_V_read => x_48_V_read,
        x_17_V_read => x_49_V_read,
        x_18_V_read => x_50_V_read,
        x_19_V_read => x_51_V_read,
        x_20_V_read => x_52_V_read,
        x_21_V_read => x_53_V_read,
        x_22_V_read => x_54_V_read,
        x_23_V_read => x_55_V_read,
        x_24_V_read => x_56_V_read,
        x_25_V_read => x_57_V_read,
        x_26_V_read => x_58_V_read,
        x_27_V_read => x_59_V_read,
        x_28_V_read => x_60_V_read,
        x_29_V_read => x_61_V_read,
        x_30_V_read => x_62_V_read,
        x_31_V_read => x_63_V_read,
        ap_return => grp_reduce_17_fu_592_ap_return,
        ap_ce => grp_reduce_17_fu_592_ap_ce);





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
                ap_return_int_reg <= add_ln703_fu_660_p2;
            end if;
        end if;
    end process;
    add_ln703_fu_660_p2 <= std_logic_vector(unsigned(grp_reduce_17_fu_592_ap_return) + unsigned(grp_reduce_17_fu_524_ap_return));
        ap_block_pp0_stage0 <= not((ap_const_boolean_1 = ap_const_boolean_1));
        ap_block_pp0_stage0_11001 <= not((ap_const_boolean_1 = ap_const_boolean_1));
        ap_block_pp0_stage0_11001_ignoreCallOp67 <= not((ap_const_boolean_1 = ap_const_boolean_1));
        ap_block_pp0_stage0_11001_ignoreCallOp68 <= not((ap_const_boolean_1 = ap_const_boolean_1));
        ap_block_state1_pp0_stage0_iter0 <= not((ap_const_boolean_1 = ap_const_boolean_1));
        ap_block_state1_pp0_stage0_iter0_ignore_call65 <= not((ap_const_boolean_1 = ap_const_boolean_1));
        ap_block_state1_pp0_stage0_iter0_ignore_call66 <= not((ap_const_boolean_1 = ap_const_boolean_1));
        ap_block_state2_pp0_stage0_iter1 <= not((ap_const_boolean_1 = ap_const_boolean_1));
        ap_block_state2_pp0_stage0_iter1_ignore_call65 <= not((ap_const_boolean_1 = ap_const_boolean_1));
        ap_block_state2_pp0_stage0_iter1_ignore_call66 <= not((ap_const_boolean_1 = ap_const_boolean_1));

    ap_return_assign_proc : process(add_ln703_fu_660_p2, ap_ce_reg, ap_return_int_reg)
    begin
        if ((ap_const_logic_0 = ap_ce_reg)) then 
            ap_return <= ap_return_int_reg;
        elsif ((ap_const_logic_1 = ap_ce_reg)) then 
            ap_return <= add_ln703_fu_660_p2;
        end if; 
    end process;


    grp_reduce_17_fu_524_ap_ce_assign_proc : process(ap_block_pp0_stage0_11001_ignoreCallOp67)
    begin
        if (((ap_const_logic_1 = ap_const_logic_1) and (ap_const_boolean_0 = ap_block_pp0_stage0_11001_ignoreCallOp67))) then 
            grp_reduce_17_fu_524_ap_ce <= ap_const_logic_1;
        else 
            grp_reduce_17_fu_524_ap_ce <= ap_const_logic_0;
        end if; 
    end process;


    grp_reduce_17_fu_592_ap_ce_assign_proc : process(ap_block_pp0_stage0_11001_ignoreCallOp68)
    begin
        if (((ap_const_logic_1 = ap_const_logic_1) and (ap_const_boolean_0 = ap_block_pp0_stage0_11001_ignoreCallOp68))) then 
            grp_reduce_17_fu_592_ap_ce <= ap_const_logic_1;
        else 
            grp_reduce_17_fu_592_ap_ce <= ap_const_logic_0;
        end if; 
    end process;

end behav;
