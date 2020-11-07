#!/usr/bin/env python
# coding: utf-8

#added time.sleep into v3
# In[13]:


from pytrends.request import TrendReq
import pandas as pd
import time

def Download(df2):
    resultX=pd.DataFrame()
    # Period-1
    dataset = []

    for x in range(0, len(df2)):
        keywords = [df2[x]]
        pytrend.build_payload(
            kw_list=keywords,
            cat=0,
            timeframe='2018-05-02 2018-12-31',
            geo='US')
        data = pytrend.interest_over_time()
        if not data.empty:
            data = data.drop(labels=['isPartial'], axis='columns')
            dataset.append(data)
        time.sleep(5)

    result1 = pd.concat(dataset, axis=1)
    #result1.to_csv('P1_search_trends_try.csv')
    resultX = pd.concat([resultX,result1])

    executionTime = (time.time() - startTime)
    print('P1_Execution time in sec.: ' + str(executionTime))

    # Period-2
    dataset = []

    for x in range(0, len(df2)):
        keywords = [df2[x]]
        pytrend.build_payload(
            kw_list=keywords,
            cat=0,
            timeframe='2019-01-01 2019-06-30',
            geo='US')
        data = pytrend.interest_over_time()
        if not data.empty:
            data = data.drop(labels=['isPartial'], axis='columns')
            dataset.append(data)
        time.sleep(5)

    result2 = pd.concat(dataset, axis=1)
    #result2.to_csv('P2_search_trends_try.csv')
    resultX = pd.concat([resultX,result2])

    executionTime = (time.time() - startTime)
    print('P2_Execution time in sec.: ' + str(executionTime))

    # Period-3
    dataset = []

    for x in range(0, len(df2)):
        keywords = [df2[x]]
        pytrend.build_payload(
            kw_list=keywords,
            cat=0,
            timeframe='2019-07-01 2019-12-31',
            geo='US')
        data = pytrend.interest_over_time()
        if not data.empty:
            data = data.drop(labels=['isPartial'], axis='columns')
            dataset.append(data)
        time.sleep(2)

    result3 = pd.concat(dataset, axis=1)
    #result3.to_csv('P3_search_trends_try.csv')
    resultX = pd.concat([resultX,result3])

    executionTime = (time.time() - startTime)
    print('P3_Execution time in sec.: ' + str(executionTime))


    resultX.to_csv('all_results_31.csv')

startTime = time.time()
pytrend = TrendReq(hl='en-US', tz=360)

colnames = ["keywords"]
df = pd.read_csv("keyword_list_31.csv", names=colnames)
df3 = df["keywords"].values.tolist()
df3.remove("Keywords")

Download(df3)









