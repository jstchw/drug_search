import './App.scss';
import AppRoutes from './AppRoutes';
import { CookieInfoFooter } from './components/CookieInfoFooter/CookieInfoFooter';
import AppProviders from './AppProviders';

const App = () => {
  return (
    <>
      <AppProviders>
        <AppRoutes />
      </AppProviders>
      <CookieInfoFooter />
    </>
  );
};

export default App;
