import React, { useState, useEffect } from 'react';
import logo from './logo.svg';
import './App.css';

function App() {

    const [currentTime, setCurrentTime] = useState(null);

    useEffect(() => {
        const interval = setInterval(() => {
            fetch('/api/time').then(res => res.json()).then(data => {
                const date = new Date(data.time * 1000);
                setCurrentTime(date.toLocaleTimeString());
            });
        }, 1000);
        return () => clearInterval(interval);
    }, []);


  return (
    <div className="App">
      <header className="App-header">
        <img src={logo} className="App-logo" alt="logo" />
        <p>The current time is {currentTime}.</p>
      </header>
    </div>
  );
}

export default App;
