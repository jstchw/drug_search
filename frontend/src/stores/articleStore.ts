import { create } from 'zustand';

interface StoreStates {
  articleTerms: {
    term: string;
    type: string;
    isRemovable: boolean;
  }[];
}

interface StoreActions {
  addArticleTerm: (term: string, type: string, isRemovable: boolean) => void;
  removeArticleTerm: (term: string, type: string) => void;
  resetArticleTerms: () => void;
}

type Store = StoreStates & StoreActions;

const useArticleStore = create<Store>((set) => ({
  articleTerms: [],
  addArticleTerm: (term: string, type: string, isRemovable: boolean) =>
    set((state) => {
      const termExists = state.articleTerms.some((item) => item.term === term);
      if (!termExists) {
        return { articleTerms: [...state.articleTerms, { term, type, isRemovable }] };
      }
      return state;
    }),
  removeArticleTerm: (term: string, type: string) =>
    set((state) => ({
      articleTerms: state.articleTerms.filter((articleTerm) => articleTerm.term !== term || articleTerm.type !== type),
    })),
  resetArticleTerms: () => set({ articleTerms: [] }),
}));

export default useArticleStore;
