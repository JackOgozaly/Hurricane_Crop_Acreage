# -*- coding: utf-8 -*-
"""
Spyder Editor

This is a temporary script file.
"""
 
import pandas as pd
print(pd.__version__)
import numpy as np
import os
import seaborn as sns

cause_of_loss_4 = pd.read_parquet('https://github.com/JackOgozaly/Hurricane_Crop_Acreage/blob/main/Data/Cause_Of_Loss/crop_loss_data_4.parquet.gzip?raw=true')

print(cause_of_loss_4)

dfd = cause_of_loss_4[(cause_of_loss_4["year_of_loss"]==2017)&(cause_of_loss_4["cause_of_loss_description"]=="Hurricane/Tropical Depression")]
dfd = dfd.drop(columns=["insurance_plan_name_abbreviation","subsidy","state/private_subsidy","additional_subsidy","efa_premium_discount","producer_paid_premium","insurance_plan_code","stage_code","net_endorsed_acres"])

