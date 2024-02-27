import { create } from 'zustand';
import { ThemeType } from 'src/types';
import Cookies from 'js-cookie';

interface StoreStates {
  theme: ThemeType;
}

interface StoreActions {
  toggleTheme: () => void;
}

type Store = StoreStates & StoreActions;

const toggleTheme = (state: Store): ThemeType => {
  const theme = state.theme === 'light' ? 'dark' : 'light';
  document.documentElement.setAttribute('data-bs-theme', theme);
  Cookies.set('theme', theme, { sameSite: 'strict' });
  return theme;
};

const useGeneralOptionsStore = create<Store>((set) => ({
  theme: (Cookies.get('theme') as ThemeType | undefined) ?? 'light',
  toggleTheme: () => set((state) => ({ theme: toggleTheme(state) })),
}));

export default useGeneralOptionsStore;
