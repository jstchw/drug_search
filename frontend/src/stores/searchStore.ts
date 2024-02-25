import { create } from 'zustand';
import { SearchOptions } from '../types';
import Cookies from 'js-cookie';
import { defaultSearchOptions } from '../constants';

interface StoreStates {
  searchOptions: SearchOptions;
  searchInput: string[];
}

interface StoreActions {
  setSearchOptions: (searchOptions: SearchOptions) => void;
  setSearchInput: (searchInput: string[]) => void;
}

type Store = StoreStates & StoreActions;

const initialSearchOptionsState = () => {
  const optionsString = Cookies.get('searchOptions');
  if (optionsString) {
    try {
      return JSON.parse(optionsString);
    } catch (e) {
      console.error(e);
    }
  }
  return defaultSearchOptions;
};

const useSearchStore = create<Store>((set) => ({
  searchOptions: initialSearchOptionsState(),
  setSearchOptions: (searchOptions: SearchOptions) => set({ searchOptions }),
  searchInput: [],
  setSearchInput: (searchInput: string[]) => set({ searchInput }),
}));

export default useSearchStore;
