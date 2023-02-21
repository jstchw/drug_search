import React, { useState, useEffect } from 'react';
import Routes from './Routes';

function App() {
    const [currentTime, setCurrentTime] = useState(null);

    useEffect(() => {
        fetch('/api/time').then(res => res.json()).then(data => {
            const date = new Date(data.time * 1000);
            setCurrentTime(date.toLocaleTimeString());
        });
    }, []);


  return (
    <div className="App">
      <header className="App-header">
        <p>The current time is {currentTime}.</p>
      </header>
    </div>
  );
}

export default App;
