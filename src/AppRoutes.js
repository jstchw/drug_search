import { BrowserRouter, Route, Routes } from "react-router-dom";
import Dashboard from "./pages/Dashboard/Dashboard";
import Search from "./pages/Search/Search";
import Account from "./pages/Account/Account";
import Settings from "./pages/Settings/Settings";

const AppRoutes = () => {
    return (
        <BrowserRouter>
            <Routes>
                <Route exact path="/" element={<Search/>} />
                <Route path="/dashboard" element={<Dashboard/>} />
                <Route path="/settings" element={<Settings/>} />
            </Routes>
        </BrowserRouter>
    );
};

export default AppRoutes;
