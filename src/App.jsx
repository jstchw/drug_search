import React, {useEffect} from 'react'
import AppRoutes from './AppRoutes'
import Cookies from 'js-cookie'
import {ThemeContext} from "./contexts/ThemeContext";

const App = () => {
    const [theme, setTheme] = React.useState(() => {
        const theme = Cookies.get('theme')
        return theme ? theme : 'light'
    })

    // whenever theme changes, update the document's theme attribute and a cookie
    useEffect(() => {
        document.documentElement.setAttribute('data-bs-theme', theme)
        Cookies.set('theme', theme)
    }, [theme])

    const toggleTheme = () => {
        setTheme(prevState => prevState === 'light' ? 'dark' : 'light')
    }

    return (
        <ThemeContext.Provider value={{theme, toggleTheme}}>
            <AppRoutes />
        </ThemeContext.Provider>
    )
}

export default App;
