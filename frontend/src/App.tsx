import './App.scss';
import AppRoutes from './AppRoutes';
import { CookieInfoFooter } from './components/CookieInfoFooter/CookieInfoFooter';
import AppProviders from './AppProviders';
import useGeneralOptionsStore from './stores/generalOptionsStore';

const initSettings = () => {
  const theme = useGeneralOptionsStore((state) => state.theme);
  document.documentElement.setAttribute('data-bs-theme', theme);
};

const App = () => {
  initSettings();
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
