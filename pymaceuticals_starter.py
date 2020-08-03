#!/usr/bin/env python
# coding: utf-8

# ## Observations and Insights 

# 

# In[25]:


# Dependencies and Setup
import matplotlib.pyplot as plt
import pandas as pd
import scipy.stats as st
import numpy as np

# files
mouse_metadata = "data/Mouse_metadata.csv"
study_results = "data/Study_results.csv"

# files to DF
mouse_df = pd.read_csv(mouse_metadata)
study_df = pd.read_csv(study_results)

#combine data frames
full_df = pd.merge(mouse_df, study_df, on = "Mouse ID", how = "outer")
full_df.head()

#sort
full_df.sort_values("Timepoint")


# In[34]:


# Checking the number of mice.
full_df['Mouse ID'].nunique()


# In[35]:


# Getting the duplicate mice by ID number that shows up for Mouse ID and Timepoint. 
find_dup = full_df[full_df.duplicated(subset=['Mouse ID','Timepoint'])]
print(find_dup)


# In[36]:


# Optional: Get all the data for the duplicate mouse ID. 
full_df.loc[full_df['Mouse ID']=='g989']


# In[37]:


# Create a clean DataFrame by dropping the duplicate mouse by its ID.
full_df = full_df[full_df['Mouse ID'] != 'g989']
full_df.head()


# In[38]:


# Checking the number of mice in the clean DataFrame.
full_df['Mouse ID'].nunique()


# ## Summary Statistics

# In[68]:


# Generate a summary statistics table of mean, median, variance, standard deviation, and SEM of the tumor volume for each regimen
full_df_sorted = full_df.sort_values(['Tumor Volume (mm3)'], ascending = True)
full_df_sorted.head()


# In[43]:


#Stats by drug grouped
drug_grouped = full_df_sorted.groupby(["Drug Regimen"])
drug_grouped

#create a variable to capture the total tumor volume for each regimen
#tumor_volume = regimen_grouped["Tumor Volume (mm3)"].sum()

#create computation for the mean of each regimen
drug_mean = drug_grouped["Tumor Volume (mm3)"].mean()

#Create computation for the median of each regimen
drug_median = drug_grouped["Tumor Volume (mm3)"].median()

#Create computation for the variance of each regimen
drug_variance = drug_grouped["Tumor Volume (mm3)"].var()

#create computation for the standard deviation of each regimen
drug_std = drug_grouped["Tumor Volume (mm3)"].std()

#create computation for the SEM
drug_sem = drug_grouped["Tumor Volume (mm3)"].sem()




# Generate a summary statistics table of mean, median, variance, standard deviation, and SEM of the tumor volume for each regimen
summary_df = pd.DataFrame({"Mean": drug_mean, "Median": drug_median, "Variance":drug_variance, "Standard Deviation": drug_std, "SEM": drug_sem})

summary_df


# In[45]:


# Use count to make datapoints
drug_datapoints = full_df.groupby(["Drug Regimen"]).count()["Mouse ID"]
drug_datapoints


# In[65]:


# generate bar plot using pyplot
drug_datapoints.plot(kind="bar", figsize=(10,5))

#set chart labels
plt.title("Data Points Visualized")
plt.xlabel("Drug Regimen")
plt.ylabel("Number of Data Points")

bar_df = full_df.sort_values('Drug Regimen')


#show chart
plt.show()
plt.tight_layout()


# In[67]:


# Generate a pie plot showing the distribution of female versus male mice using pandas

colors = ['red', 'blue']
explode = (0.1,0)
plot = mouse_gender.plot.pie(y='Total Count',figsize=(5,5), colors = colors, startangle=140, explode = explode, shadow = True, autopct="%1.1f%%")


# In[12]:


# Generate a pie plot showing the distribution of female versus male mice using pyplot


# ## Quartiles, Outliers and Boxplots

# In[69]:


# Calculate the final tumor volume of each mouse across four of the treatment regimens:  
# Capomulin, Ramicane, Infubinol, and Ceftamin
best_regimens = full_df[full_df["Drug Regimen"].isin(["Capomulin", "Ramicane", "Infubinol", "Ceftamin"])]
best_regimens = best_regimens.sort_values(["Timepoint"], ascending=True)
best_regimens

best_regimens_df = best_regimens[["Drug Regimen", "Mouse ID", "Timepoint", "Tumor Volume (mm3)"]]

best_regimens_df


# In[70]:


#Group data by Drug Regimen and Mouse ID to capture Last Tumor Measurement
best_regimens_sort = best_regimens_df.groupby(['Drug Regimen', 'Mouse ID']).last()['Tumor Volume (mm3)']
best_regimens_sort.head()

# Turn retrieved data into dataframe to easily manipulate
best_regimens_sort_df = best_regimens_sort.to_frame()
best_regimens_sort_df

#Create a list to use as labels and dataframe
top_4 = ['Capomulin', 'Ramicane', 'Infubinol','Ceftamin']

# Generate a box plot of the final tumor volume of each mouse across four regimens of interest
final_df = best_regimens_sort_df.reset_index()
tumor_lists = final_df.groupby('Drug Regimen')['Tumor Volume (mm3)'].apply(list)
tumor_list_df = pd.DataFrame(tumor_lists)
tumor_list_df = tumor_list_df.reindex(top_4)
tumor_vols = [vol for vol in tumor_list_df['Tumor Volume (mm3)']]
plt.boxplot(tumor_vols, labels=top_4)
plt.ylim(10, 80)
plt.show()
    


# In[72]:



# Generate a line plot of time point versus tumor volume for a mouse treated with Capomulin
time_vs_tumor = full_df[full_df["Mouse ID"].isin(["j119"])]
time_vs_tumor

time_vs_tumor_data = time_vs_tumor[["Mouse ID", "Timepoint", "Tumor Volume (mm3)"]]
time_vs_tumor_data

line_plot_df = time_vs_tumor_data.reset_index()
line_plot_df

line_plot_final = line_plot_df[["Mouse ID", "Timepoint", "Tumor Volume (mm3)"]]
line_plot_final

lines = line_plot_final.plot.line()


# ## Line and Scatter Plots

# In[74]:


# Generate a scatter plot of mouse weight versus average tumor volume for the Capomulin regimen
capomulin_scatter = full_df[full_df["Drug Regimen"].isin(["Capomulin"])]

capomulin_scatter_df = best_regimens[["Mouse ID","Weight (g)", "Tumor Volume (mm3)"]]

capomulin_sorted = capomulin_scatter_df.sort_values(["Weight (g)"], ascending=True)

capomulin_scatter_plot = capomulin_scatter.reset_index()

capomulin_grouped_weight = capomulin_scatter_plot.groupby("Weight (g)")["Tumor Volume (mm3)"].mean()

capo_grouped_plot = capomulin_grouped_weight.reset_index()


#capomulin_scatter = capomulin_grouped_weight.plot.scatter(x='Weight (g)', y='Tumor Volume (mm3)')
#
capomulin_scatter = capo_grouped_plot.plot(kind='scatter', x='Weight (g)', y='Tumor Volume (mm3)', grid = True, figsize= (8,8))


# ## Correlation and Regression

# In[75]:


# Calculate the correlation coefficient and linear regression model 
x_values = capo_grouped_plot["Weight (g)"]
y_values = capo_grouped_plot["Tumor Volume (mm3)"]
(slope, intercept, rvalue, pvalue, stderr) = st.linregress(x_values, y_values)

regress_values = x_values * slope + intercept
line_eq = "y =" + str(round(slope,2)) + "x + " + str(round(intercept,2))

plt.scatter(x_values, y_values)
plt.plot(x_values,regress_values,"r-")
plt.annotate(line_eq,(6,10),fontsize=10,color="red")
plt.xlabel("Weight")
plt.ylabel("Tumor Volume")
plt.title("Weight Vs. Avg Tumor Vol")
plt.show()


# for mouse weight and average tumor volume for the Capomulin regimen


# In[ ]:




