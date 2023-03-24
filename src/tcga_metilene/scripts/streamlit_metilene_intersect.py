import yfinance as yf
import streamlit as st
import pandas as pd

st.write("""
# test app

something whats happening here


""")

tickerSymbol = 'GOOGL'

tickerData = yf.Ticker(tickerSymbol)

tickerDF = tickerData.history(period='1d', start='2010-5-31', end='2020-5-31')
st.write(tickerDF)

st.line_chart(tickerDF.Close)
st.line_chart(tickerDF.Volume)
st.button('test_button')
