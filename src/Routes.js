import { Route } from "react-router-dom";
import Dashboard from "./pages/Dashboard/Dashboard";

const Routes = () => {
    return (
        <Routes>
            <Route path="/dashboard" element={Dashboard} />
        </Routes>
    );
};

export default Routes;
