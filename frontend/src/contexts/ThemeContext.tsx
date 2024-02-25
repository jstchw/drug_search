import React from 'react';
import Cookies from 'js-cookie';
import { ThemeType } from 'src/types';

interface ThemeContextProps {
  theme: ThemeType;
  toggleTheme: () => void;
}

export const ThemeContext = React.createContext<ThemeContextProps>({
  theme: 'light',
  toggleTheme: () => {},
});

export const ThemeProvider: React.FC<{ children: React.ReactNode }> = ({ children }) => {
  const [theme, setTheme] = React.useState<ThemeType>(() => {
    const theme = Cookies.get('theme') as ThemeType | undefined;
    return theme ?? 'light';
  });

  const toggleTheme = React.useCallback(() => {
    setTheme((prevState) => (prevState === 'light' ? 'dark' : 'light'));
  }, []);

  React.useEffect(() => {
    document.documentElement.setAttribute('data-bs-theme', theme);
    Cookies.set('theme', theme, { sameSite: 'strict' });
  }, [theme]);

  return <ThemeContext.Provider value={{ theme, toggleTheme }}>{children}</ThemeContext.Provider>;
};
