import React from 'react';
import { QueryClient, QueryClientProvider } from 'react-query';
import { SearchHistoryProvider } from './contexts/SearchHistoryContext';

const queryClient = new QueryClient();

interface AppProvidersProps {
  children: React.ReactNode;
}

const AppProviders: React.FC<AppProvidersProps> = ({ children }) => {
  return (
    <QueryClientProvider client={queryClient}>
      <SearchHistoryProvider>{children}</SearchHistoryProvider>
    </QueryClientProvider>
  );
};

export default AppProviders;
