import React from 'react';
import { QueryClient, QueryClientProvider } from 'react-query';
import { SearchHistoryProvider } from './contexts/SearchHistoryContext';
import { ThemeProvider } from './contexts/ThemeContext';

const queryClient = new QueryClient();

interface AppProvidersProps {
  children: React.ReactNode;
}

const AppProviders: React.FC<AppProvidersProps> = ({ children }) => {
  return (
    <QueryClientProvider client={queryClient}>
      <ThemeProvider>
        <SearchHistoryProvider>{children}</SearchHistoryProvider>
      </ThemeProvider>
    </QueryClientProvider>
  );
};

export default AppProviders;
